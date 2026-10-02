## FP64 reconstruction of FP32 cached rows, checked against an 18-real
## reference built with host complex arithmetic. A point source isolates
## individual link entries so FP32 stencil rounding cannot mask a mismatch.
import testutils
import qex, physics/qcdTypes, physics/stagGpu, backend/accel
import std/[complex, math, strformat]

qexInit()
const V = static(VLEN)
let lo = latticeFromLocalLattice(@[4,4,4,4], nRanks).newLayout
let n = lo.nSites
var rng = lo.newRngField(MRG32k3a, 987654321'u)
var g = lo.newGauge()

proc expand(s: StagGpu[V,float32]): auto =
  result = lo.newGauge()
  var u = newSeq[float32](4*s.nl*n)
  var nb = newSeq[int32](8*n)
  gpuMemCpyToCpu(addr u[0],s.lf,u.len*sizeof(float32))
  gpuMemCpyToCpu(addr nb[0],s.nbr,nb.len*sizeof(int32))
  for mu in 0..3:
    var h = newSeq[float](18*n)
    for k in 0..<n:
      let o = s.nl*mu*n+(k div V)*s.nl*V+k mod V
      let d = (k div V)*18*V+k mod V
      for c in 0..<(if s.nl==18: 18 else: 12): h[d+c*V] = float(u[o+c*V])
      if s.nl < 18:
        var rows: array[2,array[3,Complex64]]
        for a in 0..1:
          for b in 0..2: rows[a][b] = complex64(float(u[o+(6*a+2*b)*V]),float(u[o+(6*a+2*b+1)*V]))
        let phase = if s.nl==12: complex64(if nb[mu*n+k]<0: -1.0 else: 1.0)
                    else: complex64(float(u[o+12*V]),float(u[o+13*V]))
        for b in 0..2:
          let c = (b+1) mod 3
          let e = (b+2) mod 3
          let z = phase*conjugate(rows[0][c]*rows[1][e]-rows[0][e]*rows[1][c])
          h[d+(12+2*b)*V] = float(float32(z.re))
          h[d+(13+2*b)*V] = float(float32(z.im))
    copyMem(addr result[mu][0],addr h[0],h.len*sizeof(float))

proc apply[T](s: StagGpu[V,T]): seq[T] =
  var x = newSeq[T](6*n)
  if myRank==0: x[4*V] = T(1) # color 2 at even site zero
  gpuMemCpyToGpu(s.vec[5],addr x[0],x.len*sizeof(T))
  s.applyMfull(s.vec[0],s.vec[5],0.0)
  result = newSeq[T](6*n)
  gpuMemCpyToCpu(addr result[0],s.vec[0],result.len*sizeof(T))

proc difference[T](a,b: seq[T]): float =
  for i in 0..<a.len: result = max(result,abs(float(a[i])-float(b[i])))
  getDefaultComm().max(result)

suite "staggered reconstruction precision":
  for nl in [12,14,18]:
    g.random rng
    threads:
      g.projectSU
      if nl==14:
        for mu in 0..3:
          let z = newComplex(cos(0.17*float(mu+1)),sin(0.17*float(mu+1)))
          for i in g[mu]: g[mu][i] := z*g[mu][i]
      threadBarrier()
      g.setBC
      threadBarrier()
      g.stagPhase
    for fw in [0,1]:
      test &"{nl} reals, forward {fw}":
        var s = newStagGpu(g,float32,nl,fw,recon64=true)
        check s.nl==nl
        let h = expand(s)
        var r = newStagGpu(h,float32,18,fw)
        let want = r.apply()
        let wide = s.apply()
        s.recon64 = false
        let narrow = s.apply()
        let err = difference(wide,want)
        let change = difference(narrow,wide)
        echo &"reconstruction {nl} fwd {fw} reference error {err:.3e} mode difference {change:.3e}"
        check err==0
        if nl<18: check change>0
        else: check change==0
        s.recon64 = true
        check difference(s.apply(),wide)==0
        s.free()
        r.free()
        var d = newStagGpu(g,float64,nl,fw)
        let a = d.apply()
        d.recon64 = true
        check difference(a,d.apply())==0
        d.free()
qexFinalize()
