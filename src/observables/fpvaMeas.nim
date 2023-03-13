import qex
import physics/qcdTypes
import gauge
import physics/stagSolve

proc pointSource(r:Field; c:openArray[int]; ic:int) =
  let (ptRank,ptIndex) = r.l.rankIndex(c)
  threads:
    r := 0
    echo "point: ", r.norm2
    if myRank==ptRank and threadNum==0:
      r{ptIndex}[ic] := 1
    threadBarrier()
    echo "point: ", r.norm2

proc symShift(r:Field; g:Field2; x:Field; dir:int) =
  var sf = createShiftB(r, dir, 1)
  var sb = createShiftB(r, dir,-1)
  threads:
    echo "symShift: ", x.norm2
    startSB(sf, x[ix])
    startSB(sb, g[ix].adj*x[ix])
    for ir in r:
      localSB(sf, ir, r[ir] := g[ir]*it, x[ix])
      localSB(sb, ir, r[ir] += it, g[ix].adj*x[ix])
    boundarySB(sf, r[ir] += g[ir]*it)
    boundarySB(sb, r[ir] += it)
    echo "symShift: ", r.norm2

proc stagLocalMesons(v1,v2:auto, t0=0):auto =
  let l = v1.l
  let nt = l.physGeom[3]
  var c = newSeq[array[8,float]](nt)
  when true:
  #when false:
    var x1: VectorArray[3,DComplex]
    var x2: VectorArray[3,DComplex]
    for i in 0..<l.nSites:
      let t = l.coords[3][i]
      let s = (l.coords[0][i].int and 1) + ((l.coords[1][i].int and 1) shl 1) +
              ((l.coords[2][i].int and 1) shl 2)
      let tt = (t+nt-t0) mod nt
      #c[tt][s] += redot(v1{i}, v2{i})
      assign(x1, v1{i}[])
      assign(x2, v2{i}[])
      c[tt][s] += redot(x1, x2)
  else:
    threads:
      var x:VectorArray[3,SComplex]
      for i in 0..<l.nSites:
        let t = l.coords[3][i]
        let s = (l.coords[0][i].int and 1)+((l.coords[1][i].int and 1) shl 1) +
              ((l.coords[2][i].int and 1) shl 2)
        let tpar = (8*t+s) mod numThreads
        if tpar==threadNum:
          let tt = (t+nt-t0) mod nt
          #c[t][s] += v{i}.norm2()
          assign(x, v{i})
          c[tt][s] += x.norm2
  rankSum(c)
  result = c

proc sft(c:var auto, b:int) =
  for t in 0..<c.len:
    for s in 0..<8:
      if (s and b)==0:
        let c0 = c[t][s]
        let c1 = c[t][s+b]
        c[t][s] = c0 + c1
        c[t][s+b] = c0 - c1

proc printLocalMesons(c:var auto, f=1.0) =
  sft(c, 1)
  sft(c, 2)
  sft(c, 4)
  for s in 0..<8:
    echo "corner: ", s
    for t in 0..<c.len:
      let r = c[t][s]
      echo t, " ", f*r

proc `+=`[T](r:var openArray[T]; x:openArray[T]) =
  for i in 0..<r.len: r[i] += x[i]

template mysolve(dest, src) =
  #threads: dest := 0
  s.solve(dest, src, m, sp)
  threads:
    echo "dest: ", dest.norm2
    echo "dest.even: ", dest.even.norm2
    echo "dest.odd: ", dest.odd.norm2
    s.D(r, dest, m)
    threadBarrier()
    r := src - r
    threadBarrier()
    echo r.norm2

when isMainModule:
  qexInit()
  let defaultGaugeFile = "l88.scidac"
  #let defaultLat = [8,8,8,16]
  #let defaultLat = [8,8,8,8]
  #let defaultLat = [4,4,4,4]
  defaultSetup()
  var src = lo.ColorVector()
  var dest = lo.ColorVector()
  var r = lo.ColorVector()
  var destS = lo.ColorVector()
  threads:
    g.setBC
    g.stagPhase
  var s = newStag(g)
  var m = 0.1
  var sp = initSolverParams()
  sp.r2req = 1e-16
  let nt = lo.physGeom[3]
  var cl = newSeq[array[8,float]](nt)
  var cx = newSeq[array[8,float]](nt)
  var cy = newSeq[array[8,float]](nt)
  var cz = newSeq[array[8,float]](nt)
  var cs = [cx.addr,cy.addr,cz.addr]

  var t0 = 2
  let pt = [0,0,0,t0]

  for ic in 0..2:
    src.pointSource(pt, ic)
    mysolve(dest, src)
    cl += stagLocalMesons(dest, dest, t0)
    for mu in 0..2:
      symShift(r, g[mu], src, mu)
      mysolve(destS, r)
      symShift(r, g[mu], destS, mu)
      cs[mu][] += stagLocalMesons(dest, r, t0)

  let f = nt.float/(lo.physVol.float)
  printLocalMesons(cl, f)
  for mu in 0..2:
    printLocalMesons(cs[mu][], f)

  qexFinalize()
