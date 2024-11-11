import os, times

import base
import layout
import field
import qcdTypes
#import stdUtils
import solvers/cg
export cg
#import types
#import profile
#import metaUtils
import gauge/gaugeUtils

type StaggeredD*[T] = object
  sf*:seq[ShiftB[T]]
  sb*:seq[ShiftB[T]]
  sub*:string
  subset*:Subset
type Staggered*[G,T] = object
  se*,so*:StaggeredD[T]
  g*:seq[G]

template initStagDT*(l:var Layout; T:typedesc; ss:string):untyped =
  var sd:StaggeredD[T]
  sd.sf.newSeq(4)
  sd.sb.newSeq(4)
  for mu in 0..<4:
    initShiftB(sd.sf[mu], l, T, mu, 1, ss)
    initShiftB(sd.sb[mu], l, T, mu,-1, ss)
  sd.sub = ss
  sd.subset.layoutSubset(l, ss)
  sd

proc initStagD*(x:Field; sub:string):auto =
  result = initStagDT(x.l, evalType(x[0]), sub)

template initStagD3T*(l:var Layout; T:typedesc; ss:string):untyped =
  var sd:StaggeredD[T]
  sd.sf.newSeq(8)
  sd.sb.newSeq(8)
  for mu in 0..<4:
    initShiftB(sd.sf[2*mu  ], l, T, mu, 1, ss)
    initShiftB(sd.sf[2*mu+1], l, T, mu, 3, ss)
    initShiftB(sd.sb[2*mu  ], l, T, mu,-1, ss)
    initShiftB(sd.sb[2*mu+1], l, T, mu,-3, ss)
  sd.sub = ss
  sd.subset.layoutSubset(l, ss)
  sd

proc initStagD3*(x:Field; sub:string):auto =
  result = initStagD3T(x.l, evalType(x[0]), sub)

proc subdirs*(s: Staggered, dirs: seq[int]): Staggered =
  result.se.sf.newSeq(0)
  result.se.sb.newSeq(0)
  result.se.sub = s.se.sub
  result.se.subset = s.se.subset
  result.so.sb.newSeq(0)
  result.so.sf.newSeq(0)
  result.so.sub = s.so.sub
  result.so.subset = s.so.subset
  result.g.newSeq(0)
  for i in 0..<dirs.len:
    let d = dirs[i]
    result.se.sf.add s.se.sf[d]
    result.se.sb.add s.se.sb[d]
    result.so.sf.add s.so.sf[d]
    result.so.sb.add s.so.sb[d]
    result.g.add s.g[d]

proc rephase*(s: Staggered) =
  s.g.setBC
  threadBarrier()
  s.g.stagPhase

proc rephase*(s: Staggered, g: auto) =
  g.setBC
  threadBarrier()
  g.stagPhase

proc trace*(stag: Staggered, mass: float): float =
  let v = stag.g[0][0].ncols * stag.g[0].l.physVol
  result = mass*v.float

proc norm2*(stag: Staggered, mass: float): float =
  var s = 0.0
  for i in 0..<stag.g.len:
    s += stag.g[i].norm2
  let v = stag.g[0][0].ncols * stag.g[0].l.physVol
  result = mass*mass*v.float + 0.5*s


# stagDdirs
# sd, r, g, x, scales, expFlops, exp
# r = scales[mu]D_mu 

template stagDPN*(sd:openArray[StaggeredD]; r:openArray[Field];
                  g:openArray[Field2]; x:openArray[Field3];
                  expFlops:int; exp:untyped) {.dirty.} =
  tic()
  let n = sd.len
  #[
  for mu in 0..<g.len:
    for i in 0..<n:
      startSB(sd[i].sf[mu], x[i][ix])
  toc("startShiftF")
  for mu in 0..<g.len:
    for i in 0..<n:
      startSB(sd[i].sb[mu], g[mu][ix].adj*x[i][ix])
  toc("startShiftB")
  ]#
  #var rir = newAlignedMem[type(load1(r[0][0]))](n)
  #var rir:array[10,type(load1(r[0][0]))]
  #for ir{.inject.} in r[0][sd[0].subset]:
  #for ir in r[0][sd[0].subset]:
  #  for i in 0..<n:
  let ns0 = sd[0].subset.lowOuter
  let ns1 = sd[0].subset.highOuter
  let ns = ns1 - ns0
  #tFor ir, ns0, ns.pred1:
  #  for i in 0..<n:
  let nsn = ns*n
  let nr = 8
  let ni = 16
  let nin = ni * n
  let ninr = ni * nr
  let n4 = n div nr
  let n4r = n mod nr
  tFor iri, 0, nsn.pred:
      let lr = iri mod nin
      let lr64 = lr div ninr
      let lrr = lr - ninr*lr64
      var nn = nr
      if lr64>=n4: nn = n4r
      let ir = ns0 + (iri div nin)*ni + (lrr div nn)
      let i = lr64*nr + (lrr mod nn)
      #let ir = ns0 + (iri div n)
      #let i = iri mod n
      var rir{.inject,noInit.}:evalType(load1(r[i][ir]))
      exp
      for mu in 0..<g.len:
        #if mu<g.len-1:
        #  prefetchSB(sd[i].sf[mu+1], ir, x[i][ix])
        #else:
        #  prefetchSB(sd[i].sb[0], ir, x[i][ix])
        localSB(sd[i].sf[mu], ir, imadd(rir, g[mu][ir], it),load1(x[i][ix]))
      for mu in 0..<g.len:
        #if mu<g.len-1:
        #  prefetchSB(sd[i].sb[mu+1], ir, x[i][ix])
        #else:
        #  prefetchSB(sd[i].sf[mu], ir, x[i][ix])
        localSB(sd[i].sb[mu], ir, isub(rir, it), g[mu][ix].adj*x[i][ix])
      assign(r[i][ir], rir)
  toc("local", flops=n*(expFlops+g.len*(72+66+6))*sd[0].subset.len)
  #[
  for mu in 0..<g.len:
    for i in 0..<n:
      boundarySB(sd[i].sf[mu], imadd(r[i][ir], g[mu][ir], it))
  toc("boundaryF")
  for mu in 0..<g.len:
    for i in 0..<n:
      boundarySB(sd[i].sb[mu], isub(r[i][ir], it))
  toc("boundaryB")
  ]#
template stagDMN*(sd:openArray[StaggeredD]; r:openArray[Field];
                  g:openArray[Field2]; x:openArray[Field3];
                  expFlops:int; exp:untyped) =
  tic()
  let n = sd.len
  for mu in 0..<g.len:
    for i in 0..<n:
      startSB(sd[i].sf[mu], x[i][ix])
  toc("startShiftF")
  for mu in 0..<g.len:
    for i in 0..<n:
      startSB(sd[i].sb[mu], g[mu][ix].adj*x[i][ix])
  toc("startShiftB")
  var rir{.inject,noInit.}:seq[evalType(load1(r[0][0]))]
  rir.newSeq(n)
  for ir{.inject.} in r[sd.subset]:
    exp
    for mu in 0..<g.len:
      for i in 0..<n:
        localSB(sd[i].sf[mu], ir, imsub(rir[i], g[mu][ir], it),load1(x[i][ix]))
      for i in 0..<n:
        localSB(sd[i].sb[mu], ir, iadd(rir[i], it), g[mu][ix].adj*x[i][ix])
    for i in 0..<n:
      assign(r[i][ir], rir[i])
  toc("local", flops=(expFlops+g.len*(72+66+6))*sd.subset.len)
  for mu in 0..<g.len:
    for i in 0..<n:
      boundarySB(sd[i].sf[mu], imsub(r[i][ir], g[mu][ir], it))
  toc("boundaryF")
  for mu in 0..<g.len:
    for i in 0..<n:
      boundarySB(sd[i].sb[mu], iadd(r[i][ir], it))
  toc("boundaryB")

template stagDP*(sd:StaggeredD; r:Field; g:openArray[Field2];
                 x:Field3; expFlops:int; exp:untyped) =
  tic("stagDP")
  for mu in 0..<g.len:
    XoptimizeAst:
      startSB(sd.sf[mu], x[ix])
  toc("startShiftF")
  for mu in 0..<g.len:
    XoptimizeAst:
      startSB(sd.sb[mu], g[mu][ix].adj*x[ix])
  toc("startShiftB")
  #optimizeAst:
  block:
    for ir in r[sd.subset]:
      var rir{.inject,noInit.}:evalType(load1(r[ir]))
      exp
      for mu in 0..<g.len:
        #localSB(sd.sf[mu], ir, imadd(rir, g[mu][ir], it), load1(x[ix]))
        localSB(sd.sf[mu], ir, imadd(rir, g[mu][ir], it), x[ix])
      for mu in 0..<g.len:
        localSB(sd.sb[mu], ir, isub(rir, it), g[mu][ix].adj*x[ix])
      assign(r[ir], rir)
  toc("local", flops=(expFlops+g.len*(72+66+6))*sd.subset.len)
  for mu in 0..<g.len:
    template f(ir,it: untyped): untyped =
      imadd(r[ir], g[mu][ir], it)
    XoptimizeAst:
      #boundarySB(sd.sf[mu], imadd(r[ir], g[mu][ir], it))
      boundarySB2(sd.sf[mu], f)
  toc("boundaryF")
  for mu in 0..<g.len:
    template f(ir,it: untyped): untyped =
      isub(r[ir], it)
    XoptimizeAst:
      #boundarySB(sd.sb[mu], isub(r[ir], it))
      boundarySB2(sd.sb[mu], f)
  #threadBarrier()
  toc("boundaryB")

template nVecs(x:untyped):untyped =
  when compiles(nrows(x)): nrows(x)
  else: 1
template getVec(x,i:untyped):untyped = row(x,i)
template setVec(r,x,i:untyped):untyped = setRow(r,x,i)
template stagDP2*(sd:StaggeredD; r:Field; g:openArray[Field2];
                  x:Field3; expFlops:int; exp:untyped) =
  tic()
  #[
  for mu in 0..<len(g):
    startSB(sd.sf[mu], x[ix])
  toc("startShiftF")
  for mu in 0..<g.len:
    startSB(sd.sb[mu], g[mu][ix].adj*x[ix])
  toc("startShiftB")
  ]#
  const n = nVecs(x[0])
  echoImm:
    $n & ": " & $len(x[0])
  for ir{.inject.} in r[sd.subset]:
    for ic{.inject.} in 0..<n:
    #forStatic ic, 0, n.pred:
    #  block:
        var rir{.inject,noInit.}:evalType(getVec(r[ir],0))
        for mu in 0..<g.len:
          localSB(sd.sf[mu], ir, imadd(rir, g[mu][ir], it), getVec(x[ix],ic))
        for mu in 0..<g.len:
          localSB(sd.sb[mu], ir, isub(rir, it), g[mu][ix].adj*getVec(x[ix],ic))
        setVec(r[ir], rir, ic)
  toc("local", flops=n*(expFlops+g.len*(72+66+6))*sd.subset.len)
  #[
  for mu in 0..<g.len:
    boundarySB(sd.sf[mu], imadd(r[ir], g[mu][ir], it))
  toc("boundaryF")
  for mu in 0..<g.len:
    boundarySB(sd.sb[mu], isub(r[ir], it))
  #threadBarrier()
  toc("boundaryB")
  ]#
template stagDM*(sd:StaggeredD; r:Field; g:openArray[Field2];
                 x:Field3; expFlops:int; exp:untyped) =
  tic("stagDM")
  for mu in 0..<g.len:
    XoptimizeAst:
      startSB(sd.sf[mu], x[ix])
  toc("startShiftF")
  for mu in 0..<g.len:
    XoptimizeAst:
      startSB(sd.sb[mu], g[mu][ix].adj*x[ix])
  toc("startShiftB")
  for irr in r[sd.subset]:
    XoptimizeAst:
      let ir{.inject.} = irr
      var rir{.inject,noInit.}:evalType(load1(r[ir]))
      exp
      for mu in 0..<g.len:
        localSB(sd.sf[mu], ir, imsub(rir, g[mu][ir], it), load1(x[ix]))
        localSB(sd.sb[mu], ir, iadd(rir, it), g[mu][ix].adj*x[ix])
      assign(r[ir], rir)
  toc("local", flops=(expFlops+g.len*(72+66+6))*sd.subset.len)
  for mu in 0..<g.len:
    template f(ir2,it: untyped): untyped =
      imsub(r[ir2], g[mu][ir2], it)
    XoptimizeAst:
      #boundarySB(sd.sf[mu], imsub(r[ir], g[mu][ir], it))
      boundarySB2(sd.sf[mu], f)
  toc("boundaryF")
  for mu in 0..<g.len:
    template f(ir2,it: untyped): untyped =
      iadd(r[ir2], it)
    XoptimizeAst:
      #boundarySB(sd.sb[mu], iadd(r[ir], it))
      boundarySB2(sd.sb[mu], f)
  #threadBarrier()
  toc("boundaryB")

# modified D: Df + Db
template stagDFBX*(sd:StaggeredD; r:Field; g:openArray[Field2];
                   x:Field3; expFlops:int; exp:untyped) =
  tic()
  for mu in 0..<g.len:
    startSB(sd.sf[mu], x[ix])
  toc("startShiftF")
  for mu in 0..<g.len:
    startSB(sd.sb[mu], g[mu][ix].adj*x[ix])
  toc("startShiftB")
  block:
    for irr in r[sd.subset]:
      var ir{.inject.} = irr
      var rir{.inject,noInit.}:evalType(load1(r[ir]))
      exp
      for mu in 0..<g.len:
        localSB(sd.sf[mu], ir, imadd(rir, g[mu][ir], it), load1(x[ix]))
      for mu in 0..<g.len:
        localSB(sd.sb[mu], ir, iadd(rir, it), g[mu][ix].adj*x[ix])
      assign(r[ir], rir)
  toc("local", flops=(expFlops+g.len*(72+66+6))*sd.subset.len)
  for mu in 0..<g.len:
    template f(ir2,it: untyped): untyped =
      imadd(r[ir2], g[mu][ir2], it)
    boundarySB2(sd.sf[mu], f)
  toc("boundaryF")
  for mu in 0..<g.len:
    template f(ir2,it: untyped): untyped =
      iadd(r[ir2], it)
    boundarySB2(sd.sb[mu], f)
  #threadBarrier()
  toc("boundaryB")

# r = a*r + b*x + (2D)*x
proc stagD2*(sd:StaggeredD; r:SomeField; g:openArray[Field2];
             x:Field3; a:SomeNumber; b:SomeNumber2) =
  template sf0:untyped = sd.sf
  template sb0:untyped = sd.sb
  let nd = g.len
  tic()
  for mu in 0..<nd:
    optimizeAst:
      startSB(sf0[mu], x[ix])
  toc("startShiftF")
  for mu in 0..<nd:
    optimizeAst:
      startSB(sb0[mu], g[mu][ix].adj*x[ix])
  toc("startShiftB")
  for ir in r[sd.subset]:
  #let ns0 = sd.subset.lowOuter
  #let ns1 = sd.subset.highOuter
  #let ns = ns1 - ns0
  #tFor iri, 0, ns.pred:
  #  let ir = ns0 + iri
    XoptimizeAst:
      var rir{.noInit.}:evalType(r[ir])
      rir := a*r[ir] + b*x[ir]
      for mu in 0..<nd:
        localSB(sf0[mu], ir, imadd(rir, g[mu][ir], it), x[ix])
      #for mu in 0..<nd:
        localSB(sb0[mu], ir, isub(rir, it), g[mu][ix].adj*x[ix])
        #localSB(sb0[mu], ir, rir:=it, g[mu][ix].adj*x[ix])
        #var t{.noInit.}:type(load1(x[0]))
        #localSB(sb0[mu], ir, isub(rir, it), (mul(t,g[mu][ix].adj,x[ix]);t))
      assign(r[ir], rir)
  toc("local", flops=(18+nd*(72+66+6))*sd.subset.len)
  for mu in 0..<nd:
    template f(ir,it: untyped): untyped =
      imadd(r[ir], g[mu][ir], it)
    optimizeAst:
      #boundarySB(sf0[mu], imadd(r[ir], g[mu][ir], it))
      boundarySB2(sf0[mu], f)
  toc("boundaryF")
  for mu in 0..<nd:
    template f(ir,it: untyped): untyped =
      isub(r[ir], it)
    optimizeAst:
      #boundarySB(sb0[mu], isub(r[ir], it))
      boundarySB2(sb0[mu], f)
  #threadBarrier()
  toc("boundaryB")

# r = m*x + sc*D*x
proc stagDN*(sd:openArray[StaggeredD]; r:openArray[Field]; g:openArray[Field2];
             x:openArray[Field]; m:SomeNumber; sc:SomeNumber=1.0) =
  stagDPN(sd, r, g, x, 6):
    #for i in 0..<n:
    rir := m*x[i][ir]
  #r[sd.subset] := (0.5*sc)*r

# r = a*r + m*x + sc*D*x
proc stagD*(sd:StaggeredD; r:Field; g:openArray[Field2];
            x:Field; m:SomeNumber; sc:SomeNumber=1.0, a:SomeNumber=0.0) =
  stagD2(sd, r, g, x, a/(0.5*sc), m/(0.5*sc))
  r[sd.subset] := (0.5*sc)*r
  #stagDP2(sd, r, g, x, 6):
  #  #for i in 0..<n:
  #  rir := m*getVec(x[ir], ic)

proc stagD1*(sd:StaggeredD; r:Field; g:openArray[Field2];
             x:Field; m:SomeNumber) =
  stagDP(sd, r, g, x, 6):
    rir := 0

proc stagD1x*(sd:StaggeredD; r:Field; g:openArray[Field2];
              x:Field; m:SomeNumber) =
  stagDM(sd, r, g, x, 6):
    rir := 0

# r = m*x + sc*D*x
proc stagDb*(sd:StaggeredD; r:Field; g:openArray[Field2];
             x:Field; m:SomeNumber; sc:SomeNumber=1.0) =
  stagD2(sd, r, g, x, 0, m/(0.5*sc))
  #r[sd.subset] := (0.5*sc)*r
  #stagDP2(sd, r, g, x, 6):
  #  #for i in 0..<n:
  #  rir := m*getVec(x[ir], ic)

# r = m2 - Deo * Doe
proc stagD2xx*(sdx,sdy:StaggeredD; r:Field; g:openArray[Field2];
               x:Field; m2:SomeNumber) =
  tic("stagD2xx")
  var t{.global.}:evalType(x)
  if t==nil:
    threadBarrier()
    if threadNum==0:
      t = newOneOf(x)
    threadBarrier()
  #threadBarrier()
  #stagD(sdo, t, g, x, 0.0)
  #toc("init")
  block:
    stagDP(sdy, t, g, x, 0):
      rir := 0
  toc("stagDP", flops=(g.len*(72+66+6))*sdy.subset.len)
  threadBarrier()
  #toc("barrier")
  #stagD(sde, r, g, t, 0.0)
  block:
    stagDM(sdx, r, g, t, 6):
      rir := (4.0*m2)*x[ir]
  toc("stagDM", flops=(6+g.len*(72+66+6))*sdx.subset.len)
  #threadBarrier()
  #r[sde.sub] := m2*x - r
  #for ir in r[sde.subset]:
  #  msubVSVV(r[ir], m2, x[ir], r[ir])
  #r[sde.sub] := 0.25*r

proc stagD2ee*(sde,sdo:StaggeredD; r:Field; g:openArray[Field2];
               x:Field; m2:SomeNumber) =
  stagD2xx(sde, sdo, r, g, x, m2)

proc stagD2oo*(sde,sdo:StaggeredD; r:Field; g:openArray[Field2];
               x:Field; m2:SomeNumber) =
  stagD2xx(sdo, sde, r, g, x, m2)

# r = m2 + Deo * Doe
# modified D: Df + Db
proc stagD2eeFB*(sde,sdo:StaggeredD; r:Field; g:openArray[Field2];
                 x:Field; m2:SomeNumber) =
  tic()
  var t{.global.}:evalType(x)
  if t==nil:
    threadBarrier()
    if threadNum==0:
      t = newOneOf(x)
    threadBarrier()
  toc("stagD2ee init")
  block:
    stagDFBX(sdo, t, g, x, 0):
      rir := 0
  toc("stagD2ee DP")
  threadBarrier()
  #t := -t
  threadBarrier()
  toc("stagD2ee barrier")
  block:
    stagDFBX(sde, r, g, t, 6):
      rir := (4.0*m2)*x[ir]
  toc("stagD2ee DM")

#[
# r = m2 - Deo * Doe
proc stagD2eeN*(sde,sdo:StaggeredD; r:Field; g:openArray[Field2];
                x:Field; m2:SomeNumber) =
  block:
    stagDPN(sdo, t, g, x, 0):
      rir := 0
  threadBarrier()
  block:
    stagDMN(sde, r, g, t, 6):
      rir := (4.0*m2)*x[ir]
]#

proc stagPhase*(g:openArray[Field], phases:openArray[int]) =
  let l = g[0].l
  for mu in 0..<4:
    tfor i, 0..<l.nSites:
      var s = 0
      for k in 0..<4:
        s += (phases[mu] shr k) and l.coords[k][i].int
      if (s and 1)==1:
        g[mu]{i} *= -1
        #echoAll i, " ", gt[e][0,0]

template stagPhase*(g:openArray[Field]) = stagPhase(g,[8,9,11,0])

proc newStag*[G,T](g:openArray[G];v:T):auto =
  var l = g[0].l
  template t:untyped =
    evalType(v[0])
  var r:Staggered[G,t]
  r.se = initStagDT(l, t, "even")
  r.so = initStagDT(l, t, "odd")
  r.g = @g
  r

proc newStag*[G](g:openArray[G]):auto =
  var l = g[0].l
  template t:untyped =
    evalType(l.ColorVector()[0])
    #SColorVectorV
  var r:Staggered[G,t]
  r.se = initStagDT(l, t, "even")
  r.so = initStagDT(l, t, "odd")
  r.g = @g
  r

proc newStag3*[G](g:openArray[G]):auto =
  var l = g[0].l
  template t:untyped =
    evalType(l.ColorVector()[0])
  var r:Staggered[G,t]
  r.se = initStagD3T(l, t, "even")
  r.so = initStagD3T(l, t, "odd")
  r.g = @g
  r
proc newStag3*[G](g,g3:openArray[G]):auto =
  var l = g[0].l
  template t:untyped =
    evalType(l.ColorVector()[0])
  var r:Staggered[G,t]
  r.se = initStagD3T(l, t, "even")
  r.so = initStagD3T(l, t, "odd")
  var gg = newSeq[evalType(g[0])](0)
  for i in 0..<g.len:
    gg.add g[i]
    gg.add g3[i]
  r.g = gg
  r

proc D*(s:Staggered; r,x:Field; m:SomeNumber) =
  stagD(s.se, r, s.g, x, m)
  stagD(s.so, r, s.g, x, m)
proc Ddag*(s:Staggered; r,x:Field; m:SomeNumber) =
  stagD(s.se, r, s.g, x, m, -1)
  stagD(s.so, r, s.g, x, m, -1)
proc peqDdag*(s:Staggered; r,x:Field; m:SomeNumber) =
  stagD(s.se, r, s.g, x, m, -1, 1)
  stagD(s.so, r, s.g, x, m, -1, 1)
proc eoReduce*(s:Staggered; r,b:Field; m:SomeNumber) =
  # r.even = (D^+ b).even
  #dump: "b.even.norm2"
  #dump: "b.odd.norm2"
  stagD(s.se, r, s.g, b, m, -1)
  #dump: r.even.norm2
  #dump: r.odd.norm2
proc eoReconstruct*(s:Staggered; r,b:Field; m:SomeNumber) =
  # r.odd = (b.odd - Doe r.even)/m
  stagD(s.so, r, s.g, r, 0.0, -1.0/m)
  threadBarrier()
  r.odd += b/m

# (d/dg) redot[ (2D)*x, c ]
proc stagD2deriv*(s:Staggered; g:openArray[Field]; c:Field2; x:Field3) =
  template sef:untyped = s.se.sf
  template sof:untyped = s.so.sf
  let nd = g.len
  tic()
  for mu in 0..<nd:
    startSB(sef[mu], x[ix])
  for mu in 0..<nd:
    startSB(sof[mu], x[ix])
  toc("startShiftF")
  for ir in g[0][s.se.subset]:
    for mu in 0..<nd:
      var gmu = g[mu]
      localSB(sef[mu], ir, peqOuter(gmu[ir][], c[ir][], it[]), x[ix])
  for ir in g[0][s.so.subset]:
    for mu in 0..<nd:
      var gmu = g[mu]
      localSB(sof[mu], ir, peqOuter(gmu[ir][], c[ir][], it[]), x[ix])
  toc("local")
  for mu in 0..<nd:
    var gmu = g[mu]
    template f(ir,it: untyped): untyped =
      peqOuter(gmu[ir][], c[ir][], it[])
    boundarySB2(sef[mu], f)
    boundarySB2(sof[mu], f)
  toc("boundaryF")
  threadBarrier()
  for mu in 0..<nd:
    startSB(sef[mu], c[ix])
  for mu in 0..<nd:
    startSB(sof[mu], c[ix])
  toc("startShiftF")
  for ir in g[0][s.se.subset]:
    for mu in 0..<nd:
      localSB(sef[mu], ir, meqOuter(g[mu][ir], x[ir], it), c[ix])
  for ir in g[0][s.so.subset]:
    for mu in 0..<nd:
      localSB(sof[mu], ir, meqOuter(g[mu][ir], x[ir], it), c[ix])
  toc("local")
  for mu in 0..<nd:
    template f(ir,it: untyped): untyped =
      meqOuter(g[mu][ir], x[ir], it)
    boundarySB2(sef[mu], f)
    boundarySB2(sof[mu], f)
  toc("boundaryF")

proc stagDeriv*(s:Staggered; f:openArray[Field]; x:Field2) =
  template sef:untyped = s.se.sf
  template sof:untyped = s.so.sf
  let nd = f.len
  tic()
  for mu in 0..<nd:
    startSB(sef[mu], x[ix])
  for mu in 0..<nd:
    startSB(sof[mu], x[ix])
  toc("startShiftF")
  for ir in f[0][s.se.subset]:
    for mu in 0..<nd:
      var gmu = f[mu]
      localSB(sef[mu], ir, peqOuter(gmu[ir][], x[ir][], it[]), x[ix])
  for ir in f[0][s.so.subset]:
    for mu in 0..<nd:
      var gmu = f[mu]
      localSB(sof[mu], ir, meqOuter(gmu[ir][], x[ir][], it[]), x[ix])
  toc("local")
  for mu in 0..<nd:
    var gmu = f[mu]
    template fp(ir,it: untyped): untyped =
      peqOuter(gmu[ir][], x[ir][], it[])
    boundarySB2(sef[mu], fp)
    template fm(ir,it: untyped): untyped =
      meqOuter(gmu[ir][], x[ir][], it[])
    boundarySB2(sof[mu], fm)
  toc("boundaryF")
  threadBarrier()
  s.rephase f

#[
template foldl*(f,n,op:untyped):untyped =
  var r:evalType(f(0))
  r = f(0)
  for i in 1..<n:
    let
      a {.inject.} = r
      b {.inject.} = f(i)
    r = op
  r
]#

when isMainModule:
  import rng, strutils, seqUtils
  proc runtest(v1,v2,sdAll,sdEven,sdOdd,s,m:any) =
    let g = s.g
    let lo = g[0].l
    const nv = nVecs(v1[0])
    threads:
      v1 := 0
      #v2 := 1
      if myRank==0 and threadNum==0:
        when compiles(v1[0].len):
          v1{0}[0] := 1
        else:
          v1{0} := 1
      threadBarrier()
      echo v1.norm2

      stagDb(sdAll, v2, g, v1, m)
      threadBarrier()
      echo v2.norm2
      #echo v2
      s.D(v2, v1, m)
      threadBarrier()
      echo v2.norm2

      for e in v1:
        template x(d:int):untyped = lo.vcoords(d,e)
        when compiles(v1[e].len):
          v1[e][0].re := foldl(x, 4, a*10.int16+b)
        else:
          for i in 0..<v1[e].ncols:
            v1[e][0,i].re := foldl(x, 4, a*10.int16+b)
        #echo v1[e][0]
      threadBarrier()
      stagD(sdAll, v2, g, v1, 0.5)
      echo v1[0][0]
      echo v2[0][0]

    #let nrep = int(1e7/lo.physVol.float)
    let nrep = int(2e8/lo.physVol.float)
    #let nrep = int(1e9/lo.physVol.float)
    #let nrep = 1
    template makeBench(name:untyped; bar:untyped):untyped =
      proc `name T`(sd,v1,v2:any, ss="all") =
        resetTimers()
        var t0 = epochTime()
        threads:
          for rep in 1..nrep:
            stagDb(sd, v2, g, v1, 0.5)
            when bar: threadBarrier()
        var t1 = epochTime()
        let dt = t1-t0
        #var vol = lo.physVol.float
        var vol = lo.nSites.float
        if sd.sub != "all": vol *= 0.5
        let flops = nv * (6.0+g.len*2.0*72.0) * vol
        echo ss & "secs: ", dt, "  mf: ", (nrep.float*flops)/(1e6*dt)
        echoTimers()
      template name(sd:any, ss="all") = `name T`(sd, v1, v2, ss)
    subst(bench,_,benchB,_):
      makeBench(bench, false)
      makeBench(benchB, true)
      bench(sdAll, "all  ")
      benchB(sdAll, "all  ")
      bench(sdEven, "even ")
      benchB(sdEven, "even ")
      bench(sdOdd, "odd  ")
      benchB(sdOdd, "odd  ")
    proc benchEO() =
      resetTimers()
      var t0 = epochTime()
      threads:
        for rep in 1..nrep:
          stagD2ee(sdEven, sdOdd, v2, g, v1, 0.1)
      var t1 = epochTime()
      let dt = t1-t0
      #var vol = 0.5 * lo.physVol.float
      var vol = 0.5 * lo.nSites.float
      let flops = nv * (6.0+g.len*2.0*2.0*72.0) * vol
      echo "EO   secs: ", dt, "  mf: ", (nrep.float*flops)/(1e6*dt)
      #echoTimers()
    benchEO()

  qexInit()
  echo "rank ", myRank, "/", nRanks
  let cp = commandLineParams()
  #var lat = [4,4,4,4]
  var lat = [8,8,8,8]
  #var lat = [8,8,8,16]
  #var lat = [8,8,16,16]
  #var lat = [16,16,16,8]
  #var lat = [16,16,16,16]
  #var lat = [16,16,16,32]
  if cp.len>0:
    var i0 = 0
    if cp[0][0] notin {'0'..'9'}: inc i0
    for i in 0..<lat.len:
      lat[i] = (if (i0+i)<cp.len: parseInt(cp[i0+i]) else: lat[i-1])
  var lo = newLayout(lat)
  var v1 = lo.ColorVector()
  var v2 = lo.ColorVector()
  var rs = newRNGField(RngMilc6, lo, intParam("seed", 987654321).uint64)
  var g:array[4,type(lo.ColorMatrix())]
  for i in 0..<4:
    g[i] = lo.ColorMatrix()
    threads:
      g[i] := 1
  threads:
    g.setBC
    threadBarrier()
    for i in 0..<4:
      echo g[i].norm2
    threadBarrier()
    g.stagPhase
    threadBarrier()
    for i in 0..<4:
      echo g[i].norm2

  #g.loadGauge("l88.scidac")
  var sdAll = initStagD(v1, "all")
  var sdEven = initStagD(v1, "even")
  var sdOdd = initStagD(v1, "odd")
  var s = newStag(@g)
  var m = 0.1
  echo "done newStag"

  runtest(v1, v2, sdAll, sdEven, sdOdd, s, m)
  echoTimers()

  var sdAll3 = initStagD3(v1, "all")
  var sdEven3 = initStagD3(v1, "even")
  var sdOdd3 = initStagD3(v1, "odd")
  var g3:array[8,type(lo.ColorMatrix())]
  for i in 0..3:
    g3[2*i  ] = g[i]
    #g3[2*i+1] = g[i]
    g3[2*i+1] = lo.ColorMatrix()
    g3[2*i+1].randomSU rs
  var s3 = newStag3(@g3)

  runtest(v1, v2, sdAll3, sdEven3, sdOdd3, s3, m)

  #[
  const nc = v1[0].len
  const nr = 8
  type MX* = Field[VLEN,MatrixArray[nr,nc,SComplexV]]
  #type MX* = Field[VLEN,MatrixArray[nr,nc,DComplexV]]
  var m1,m2: MX
  m1.new(lo)
  m2.new(lo)
  var sdAllM = initStagD(m1, "all")
  var sdEvenM = initStagD(m1, "even")
  var sdOddM = initStagD(m1, "odd")
  var sM = newStag(@g,m1)
  echo "testing multi matrix: ", nr
  stagD(sdAllM, m2, g, m1, m)
  runtest(m1, m2, sdAllM, sdEvenM, sdOddM, sM, m)
  echoTimers()
  ]#

  #[
  #const n = 4
  var n = 4
  if cp.len>0:
    for i in 0..<cp.len:
      if cp[i][0]=='n':
        n = parseInt(cp[i][1..^1])
        break
  echo "n: ", n
  var v1a = newSeq[type(v1)](n)
  var v2a = newSeq[type(v2)](n)
  var sda = newSeq[type(sdAll)](n)
  var sda3 = newSeq[type(sdAll3)](n)
  #var sa = array[n,type(s)]
  v1a[0] = v1
  v2a[0] = v2
  sda[0] = sdAll
  sda3[0] = sdAll3
  #sda[0] = sdEven
  #sa[0] = s
  for i in 1..<n:
    v1a[i] = lo.ColorVector()
    v1a[i] := 1
    v2a[i] = lo.ColorVector()
    sda[i] = initStagD(v1, "all")
    sda3[i] = initStagD3(v1, "all")
    #sa[i] = newStag(@g)

  let nrep = int(2e7/lo.physVol.float)
  template makeBenchN(name:untyped; bar:bool):untyped =
    proc name(sd,g:any, ss="all") =
      resetTimers()
      var t0 = epochTime()
      threads:
        for rep in 1..nrep:
          stagDN(sd, v2a, g, v1a, 0.5)
          when bar: threadBarrier()
      var t1 = epochTime()
      let dt = t1-t0
      #var vol = lo.physVol.float
      var vol = lo.nSites.float
      if sd[0].sub != "all": vol *= 0.5
      let flops = n*(6.0+g.len*2.0*72.0) * vol
      echo ss & "secs: ", dt, "  mf: ", (nrep.float*flops)/(1e6*dt)
      #echoTimers()

  makeBenchN(benchN, false)

  stagDN(sda, v2a, g, v1a, m)
  benchN(sda, g)
  #echoTimers()

  #stagDN(sda3, v2a, g3, v1a, m)
  #benchN(sda3, g3)
  ]#

  qexFinalize()
