## Original version of src/observables/nlStagMesonSpec.nim
## Author: Xiao-Yong Jin

#[
  Port of nlstagspect.lua

  Compute non–local staggered meson spectroscopy correlators
  following M. Golterman (1986).  The program
  builds 8 corner wall sources on a chosen source time–slice,
  inverts to obtain propagators and then contracts them to the
  20 time–local meson operators.
]#

import qex, gauge/hypsmear, gauge/gaugefix, physics/stagSolve
import times, strformat, os, math

# XOR helper that accepts any number of integer arguments.
func ps*(a: varargs[int]): int =
  for x in a:
    result = result xor x

# Extract 0/1 corner bits.
template px*(i:int):int = (i and 1)
template py*(i:int):int = ((i shr 1) and 1)
template pz*(i:int):int = ((i shr 2) and 1)
template pt*(i:int):int = ((i shr 3) and 1)

# Convert a list of directions (1=x,2=y,3=z,4=t) into the corner index.
func shiftIdx*(idx: int; delta: seq[int]): int =
  result = idx
  for d in delta:
    result = result xor (1 shl d)

proc shiftAndRephase*(dest: Field; src: Field2; delta: seq[int]; phase: int; base: array[4,int]) =
  ## dest = shifted(src) * phase
  #echo "shift: ",delta,"  rephase: ",phase,"  rel(xyzt): ",base
  threads:
    dest := src

  for mu in delta:
    var sf = newShifter(dest, mu, 1)
    var sb = newShifter(dest, mu, -1)
    discard sf ^* dest
    discard sb ^* dest
    threads:
      # TODO: normalize
      dest := sf.field + sb.field

  let bits = [phase.px==1, phase.py==1, phase.pz==1, phase.pt==1]

  threads:
    for i in dest.sites:
      var expo = 0
      for mu in 0..3:
        if bits[mu]:
          expo += dest.l.coords[mu][i].int - base[mu]
      if (expo and 1) == 1:
        dest{i} *= -1

# inner‑product (real part) accumulated over the lattice and binned by t‑slice
proc dotByTimeslice*(v1,v2:auto): seq[float] =
  let nt = v1.l.physGeom[^1]
  result = newSeq[float](nt)
  for i in v1.sites:
    let t = v1.l.coords[3][i]
    result[t] += v1{i}.redot v2{i}
  result.rankSum

# -------------------------------------------------------------------
# Enumerate the 20 meson operators
# - Golterman, 1986
# - Ishizuka, Fukugita, Mino, Okawa, Ukawa, 1994
# IMPORTANT: We are counting from 0.
# -------------------------------------------------------------------

type MesonTerm* = tuple[delta: seq[int], phase: int, factor: int]

const eta  = [0, 1, 3, 7]
const zeta = [14,12, 8, 0]
const eps  = 15

const klm = [(0,1,2),(1,2,0),(2,0,1)]
const klma = [(0,1,2),(0,2,1),(1,0,2),
              (1,2,0),(2,0,1),(2,1,0)]

proc buildMesonOps*(): seq[seq[MesonTerm]] =
  ## Return 20 meson ops, each has one or more terms, total 64
  # TODO: label the states in the output
  result.newSeq(20)
  template add(n:int, d:seq[int], ph:int, f:int=1) =
    result[n].add (d, ph, f)

  add(0, @[], 0)
  add(1, @[], ps(eta[3], zeta[3]))
  add(16, @[0,1,2], ps(eta[0], eta[1], eta[2]))
  add(17, @[0,1,2], ps(eta[3], zeta[3], eta[0], eta[1], eta[2]))
  for k in 0..2:
    add(2, @[], ps(eta[k], zeta[k], eps))
    add(3, @[], ps(eta[3], zeta[3], eta[k], zeta[k], eps))
    add(4, @[k], eta[k])
    add(5, @[k], ps(eta[3], zeta[3], eta[k]))
    add(6, @[k], ps(zeta[k], eps))
    add(7, @[k], ps(eta[3], zeta[3], zeta[k], eps))
  for arr in klm:
    # let k = arr[0]; let l = arr[1]
    let (k,l,_) = arr
    add(10, @[k,l], ps(eta[k], eta[l]), 2)    # factor arises from epsilon_klm
    add(11, @[k,l], ps(eta[3], zeta[3], eta[k], eta[l]), 2)
    add(12, @[k,l], ps(zeta[k], zeta[l]), 2)
    add(13, @[k,l], ps(eta[3], zeta[3], zeta[k], zeta[l]), 2)
    add(18, @[0,1,2], ps(eta[k], zeta[k], eps, eta[0], eta[1], eta[2]))
    add(19, @[0,1,2], ps(eta[3], zeta[3], eta[k], zeta[k], eps, eta[0], eta[1], eta[2]))
  for arr in klma:
    # let k = arr[0]; let l = arr[1]; let m = arr[2]
    let (k,l,m) = arr
    add(8, @[l], ps(eta[k], zeta[k], eta[l], eps))
    add(9, @[l], ps(eta[3], zeta[3], eta[k], zeta[k], eta[l], eps))
    add(14, @[k,l], ps(eta[m], zeta[m], eta[k], zeta[l]))
    add(15, @[k,l], ps(eta[3], zeta[3], eta[m], zeta[m], eta[k], zeta[l]))

when isMainModule:
  qexInit()

  letParam:
    inlat = ""                     # input gauge configuration (SciDAC)
    mass  = 0.1                    # staggered quark mass
    srcT  = 0                      # source time–slice

    resid = 1e-10                  # squared residual tolerance
    cg_max = 100_000               # maximum CG iterations

    fixGauge: bool = true
    gfstop = 1e-8                  # gauge fixing stopping condition
    gforf = 1.75                   # gauge fixing over relaxation factor
    gfoutlat = ""                  # if non‑empty save gauge‑fixed config to this file

    showTimers: bool = 0           # print timers when finished

  echoParams()
  echo "rank ", myRank, "/", nRanks
  threads: echo "thread ", threadNum, "/", numThreads

  if inlat.len == 0 or not fileExists(inlat):
    qexError "Invalid gauge file, inlat = '", inlat, "'."

  let lat = getFileLattice inlat

  var lo = lat.newLayout
  var g = lo.newGauge

  qexLog "Start loading configuration."

  echo "latsize: ",lo.physGeom
  echo "volume: ",lo.physVol

  if g.loadGauge(inlat) != 0:
    qexError "Failed to load gauge field '", inlat, "'."

  qexLog "Finished loading conifguration."

  threads:
    block:
      let d = g.checkSU
      echo "unitary deviation avg: ", d.avg, "  max: ", d.max
    threadBarrier()
    g.projectSU
    threadBarrier()
    block:
      let d = g.checkSU
      echo "new unitary deviation avg: ", d.avg, "  max: ", d.max
  g.echoPlaq

  qexLog "Finished re-unitarization."

  if fixGauge:
    tic "Coulomb gauge fixing"

    var tmat = lo.ColorMatrix()
    threads: tmat := 1
    getGaugeFixTransform(tmat, g, @[0,1,2], gfstop, gforf, verb=0)

    var gFix = lo.newGauge
    gFix.gaugeTransform(g, tmat)
    g = gFix

    qexLog "Coulomb gauge fixing done. Time: " ,getElapsedTime()

    if gfoutlat.len > 0:
      echo "Saving gauge‑fixed configuration to '", gfoutlat, "' …"
      discard g.saveGauge(gfoutlat, "D")

  qexGC("after gauge fixing")

  threads:
    g.setBC
    threadBarrier()
    g.stagPhase

  let stag = g.newStag

  var sp = initSolverParams()
  sp.r2req = resid*resid
  sp.maxits = cg_max

  let mesonOps = buildMesonOps()

  var props: array[8, typeof(lo.ColorVector)]

  var nMesonTot = 0
  for ops in mesonOps.items:
    if ops.len>0: nMesonTot += ops.len

  let ntTotal = lat[^1]
  var mesons = newSeq[seq[float]](nMesonTot)
  for i in 0..<nMesonTot:
    mesons[i].newSeq(ntTotal)

  # Helper to map (operator n, local index i) → global meson index
  proc mesonIndex(n, i:int):int =
    var off = 0
    for k in 0..<n:
      off += mesonOps[k].len
    result = off + i

  for c in 0..<g[0][0].nrows:
    tic "for one color"

    for corner in 0..7:
      var src = lo.ColorVector()
      threads:
        src := 0
        threadBarrier()
        for i in src.sites:
          let cc = lo.coords
          if cc[3][i] == srcT and
            (cc[0][i] and 1) == px(corner) and
            (cc[1][i] and 1) == py(corner) and
            (cc[2][i] and 1) == pz(corner):
            src{i}[c] := 1.0
        threadBarrier()
        echo "color: ",c,"  time slice: ",srcT,"  corner: ",corner,"  src norm2: ",src.norm2

      props[corner] = lo.ColorVector()
      threads:
        props[corner] := 0
      stag.solve(props[corner], src, mass, sp)
      threads:
        echo "sol norm2: ",props[corner].norm2
    qexLog "Solves for color ",c," done. Time: ",getElapsedTime()

    # TODO: indexing is convoluted
    var midx = 0
    for n in 0..<mesonOps.len:
      let opList = mesonOps[n]
      for termIdx, term in opList:
        let globalIdx = mesonIndex(n, termIdx)

        inc midx
        #echo "color: ",c,"  op: ",n,"  meson_idx: ",midx
        for forward in 0..7:
          let backward = shiftIdx(forward, term.delta)
          #echo "combine prop: ",forward," & ",backward

          # Build shifted + rephased copy of `props[backward]`.
          var tmp = lo.ColorVector()
          let base = [px(forward), py(forward), pz(forward), srcT]
          shiftAndRephase(tmp, props[backward], term.delta,
                          15 - term.phase, base)

          # Correlator for this pair, binned by timeslice
          let corr = dotByTimeslice(props[forward], tmp)

          for t in 0..<ntTotal:
            mesons[globalIdx][t] -= term.factor.float * corr[t]
    qexLog "Contraction for color ",c," done."

  # TODO: printing is convoluted
  var idx = 0
  for n in 0..<mesonOps.len:
    for _ in 0..<mesonOps[n].len:
      echo "BEGIN MESON ", idx
      let d = mesons[idx]
      for dt in 0..<ntTotal:
        let i = (srcT + dt) mod ntTotal
        echo fmt"{dt}	{d[i]}"
      echo "END MESON ", idx
      inc idx

  if showTimers: echoTimers()
  qexFinalize()