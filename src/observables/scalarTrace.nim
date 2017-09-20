#[
  Port of scalar propagator code from FUEL.

  Only 4D SU(3) for now.

  Original comments:

  JXY 2017-09-20
  Scale trace by 1/Nc
  output PBP spatial sum normalized by SpatialVolume

  JXY 2016-12-07
  Combine even-odd and corner to one file
  Incorporating James's single file writer.

  ESW 2016-11-14
  Script file that generates noise sources, dilutes in time and corners, then saves packed traces.
  Optionally, this script saves the original noise sources and propagators.
]#

import base, layout, field

# Move dilution support to a separate module once more codes depend on it.
# Or generalize/integrate it with subset support.
type
  DilutionKind = enum
    dkEvenOdd, dkCorners3D
  Dilution = object
    case kind: DilutionKind
    of dkEvenOdd: eo: range[0..1]
    of dkCorners3D: c3d: range[0..7]
proc `$`(x: Dilution): string =
  case x.kind
  of dkEvenOdd: "EvenOdd " & $x.eo
  of dkCorners3D: "Corners3D " & $x.c3d

iterator sites(l: Layout, d: Dilution): int =
  case d.kind
  of dkEvenOdd:
    # Assuming even-odd layout
    if d.eo == 0: itemsI(0, l.nEven)
    else: itemsI(l.nEven, l.nSites)
  of dkCorners3D:
    let
      n = l.nSites
      a = (threadNum*n) div numThreads
      b = (threadNum*n+n) div numThreads
      c = d.c3d
    var i = a
    while i < b:
      if ((l.coords[0][i].int and 1) +
         ((l.coords[1][i].int and 1) shl 1) +
         ((l.coords[2][i].int and 1) shl 2)) == c:
        yield i
        i.inc

iterator dilution(dl:string): Dilution =
  case dl
  of "EO":
    yield Dilution(kind:dkEvenOdd, eo:0)
    yield Dilution(kind:dkEvenOdd, eo:1)
  of "CORNER":
    for i in 0..7:
      yield Dilution(kind:dkCorners3D, c3d:i)
  else:
    echo "ERROR: unsupported dilution type: ",dl
    qexAbort()

when isMainModule:
  import qex, gauge, physics/qcdTypes, physics/stagD, gauge/hypsmear
  import times

  qexInit()

  # Accept command line arguments and set up parameters
  let
    # Required parameters
    inlat = strParam("inlat")   # Input gauge file name
    outfn = strParam("outfn")   # Save results to outfn.{trace,noise,prop}
    mass = floatParam("mass")   # The quark mass

    # optional parameters
    cg_prec = floatParam("cg_prec", 1e-9) # Max residual with default
    cg_max = intParam("cg_max", 100_000)  # Max number of iterations
    num_stoch = intParam("num_stoch", 1)  # Number of stochastic sources
    improved_trace = intParam("improved_trace", 1).bool
    save_props = intParam("save_props", 0).bool
    source_type = strParam("source_type", "Z4") # Z4, Z2, U1, Gauss
    dilute_type = strParam("dilute_type", "EO") # EO, CORNER
    seed = intParam("seed", int(1000*epochTime())).uint64
    # write_group_size = intParam("write_group_size", 128)

  var sp = initSolverParams()
  sp.r2req = cg_prec
  sp.maxits = cg_max

  # Load lattice and determine size and parameters
  template getLat(fn:string): seq[int] =
    if inlat.len == 0:
      echo "WARNING: no input lattice file specicified."
      @[8,8,8,8]
    else:
      let lat = inlat.getFileLattice
      if lat.len == 0:
        echo "ERROR: getFileLattice failed on '", inlat, "'."
        qexAbort()
      lat
  let
    lat = inlat.getLat
    nt = lat[^1]
  var
    lo = lat.newLayout
    g = lo.newGauge
    r = newRNGField(RngMilc6, lo, seed)
  if inlat.len == 0:
    threads: g.random r
  else:
    if 0 != g.loadGauge(inlat):
      echo "ERROR: loadGauge failed on '", inlat, "'."
      qexAbort()

  echo "latsize = ",lo.physGeom
  echo "volume = ",lo.physVol
  echo "mass = ",mass
  echo "seed = ",seed
  echo "cg_prec = ",cg_prec
  echo "num_stoch = ",num_stoch
  echo "source_type = ",source_type
  echo "dilute_type = ",dilute_type

  var spatv = 1
  for i in 0..<lat.len-1: spatv *= lat[i]

  proc printPlaq(g: any) =
    let
      p = g.plaq
      sp = 2.0*(p[0]+p[1]+p[2])
      tp = 2.0*(p[3]+p[4]+p[5])
    echo "plaq ",p
    echo "plaq ss: ",sp," st: ",tp," tot: ",p.sum

  # XXX maybe check unitarity deviation?
  # This is giving NaNs
  #[
  threads:
    for mu in 0..<g.len:
      for i in g[mu]:
        g[mu][i].projectU g[mu][i]
  ]#

  g.printPlaq

  var
    info: PerfInfo
    coef = HypCoefs(alpha1:0.4, alpha2:0.5, alpha3:0.5)
  echo "smear = ",coef
  var sg = lo.newGauge
  coef.smear(g, sg, info)

  threads:
    sg.setBC
    sg.stagPhase
  var s = sg.newStag

  var
    eta = lo.ColorVector        # Noise source
    phi = lo.ColorVector        # Propagator
    trce = lo.Complex     # Trace
    tmps = lo.ColorVector

  for i in 0..<num_stoch:
    echo "Starting work on noise source ",i," ."
    phi := 0
    trce := 0

    echo "Generating a ",source_type," noise source."
    case source_type
    of "Z4": eta.z4 r
    of "Z2": eta.z2 r
    of "U1": eta.u1 r
    of "Gauss": eta.gaussian r
    else:
      echo "ERROR: Invalid noise type ",source_type,"."
      qexAbort()

    if save_props:              # XXX save eta
      echo "WARNING: save not implemented"

    for t in 0..<nt:
      for dl in dilution(dilute_type):
        echo "Source ",i," dilution pattern ",dl," timeslice ",t," ."
        tmps := 0

        # XXX implement subset later
        for i in lo.sites(dl):
          if lo.coords[^1][i] == t:
            # tmps{i} := eta{i}    # Doesn't work
            forO c, 0, tmps{0}.len-1:
              tmps{i}[c].re := eta{i}[c].re
              tmps{i}[c].im := eta{i}[c].im
        phi := 0
        s.solve(phi, tmps, mass, sp)

        if save_props:        # XXX save phi
          echo "WARNING: save not implemented"

        if improved_trace:
          echo "Computing the improved trace."
          for i in trce:
            trce[i] += mass * phi[i].dot phi[i]
        else:
          echo "Computing the unimproved trace."
          for i in trce:
            trce[i] += tmps[i].dot phi[i]

    let invnc = 1.0 / tmps[0].len.float
    trce *= invnc

    var est = newseq[float](nt)
    for i in lo.sites:
      var t:float
      t := trce{i}.re
      est[lo.coords[3][i]] += t
    est.ranksum
    for t in 0..<nt:
      echo "initsrc ",i," timeslice ",t," pbp ",est[t] / spatv.float

  # XXX save trce

  qexFinalize()
