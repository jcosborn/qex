## Independent dense references shared by numerical and graph SU(3) tests.
import qex
import std/[math, strutils]

type
  M = MatrixArray[3,3,ComplexType[float64]]
  A = MatrixArray[8,8,float64]
  V = VectorArray[8,float64]
  R[n:static int] = array[n,array[n,float64]]
  C = R[6]  # Real block representation of a complex 3x3 matrix.
  FiniteRef* = object
    order*, scale*: int
    phi*, dphi*, k*: A
    grad*, pullback*: M
    apply*: V
    log*, pair*, dlog*, invNorm*, cond*: float
  CaseRef* = object
    name*: string
    m*, dm*, c*: M
    x*, dx*, d*, b*, analyticPhi*: A
    v*: V
    fNorm*, dNorm*, analyticLog*, analyticDLog*: float
    analyticInvNorm*, analyticCond*, selectedLog*: float
    finite*: seq[FiniteRef]
  ExpRef* = object
    name*: string
    f*, dm*, alpha*, c*, e*, de*, da*: M
    branch*: int
    checkBranch*: bool
    field*, alphaPair*: float
  References* = object
    cases*: seq[CaseRef]
    exp*, simd*: seq[ExpRef]

proc eye[n:static int](): R[n] =
  for i in 0..<n: result[i][i] = 1
proc `+`[n:static int](a,b: R[n]): R[n] =
  for i in 0..<n:
    for j in 0..<n: result[i][j] = a[i][j]+b[i][j]
proc `-`[n:static int](a,b: R[n]): R[n] =
  for i in 0..<n:
    for j in 0..<n: result[i][j] = a[i][j]-b[i][j]
proc `*`[n:static int](s: float, a: R[n]): R[n] =
  for i in 0..<n:
    for j in 0..<n: result[i][j] = s*a[i][j]
proc mm[n:static int](a,b: R[n]): R[n] =
  for i in 0..<n:
    for j in 0..<n:
      var x = 0.0
      for k in 0..<n: x += a[i][k]*b[k][j]
      result[i][j] = x
proc tp[n:static int](a,b: R[n]): float =
  for i in 0..<n:
    for j in 0..<n: result += a[i][j]*b[j][i]
proc dotRef[n:static int](a,b: R[n]): float =
  for i in 0..<n:
    for j in 0..<n: result += a[i][j]*b[i][j]
proc normRef[n:static int](a: R[n]): float = sqrt(dotRef(a,a))
proc trans[n:static int](a: R[n]): R[n] =
  for i in 0..<n:
    for j in 0..<n: result[i][j] = a[j][i]
proc cm(re,im: R[3]): C =
  for i in 0..<3:
    for j in 0..<3:
      result[i][j] = re[i][j]
      result[i+3][j+3] = re[i][j]
      result[i+3][j] = im[i][j]
      result[i][j+3] = -im[i][j]
proc cnorm(a:C):float = sqrt(0.5*dotRef(a,a))
proc tah(a:C):C =
  result = 0.5*(a-trans(a))
  let t = (result[3][0]+result[4][1]+result[5][2])/3
  for i in 0..<3:
    result[i+3][i] -= t
    result[i][i+3] += t
proc toMat(a:C):M =
  for i in 0..<3:
    for j in 0..<3:
      result[i,j].re = a[i][j]
      result[i,j].im = a[i+3][j]
proc toAdj(a:R[8]):A =
  for i in 0..<8:
    for j in 0..<8: result[i,j] = a[i][j]
proc toVec(a:array[8,float64]):V =
  for i in 0..<8: result[i] = a[i]

proc generators(): array[8,C] =
  for k in 0..<8:
    var re,im:R[3]
    case k
    of 0,3,5:
      let (i,j) = if k==0: (0,1) elif k==3: (0,2) else: (1,2)
      im[i][j] = -0.5; im[j][i] = -0.5
    of 1,4,6:
      let (i,j) = if k==1: (0,1) elif k==4: (0,2) else: (1,2)
      re[i][j] = -0.5; re[j][i] = 0.5
    of 2:
      im[0][0] = -0.5; im[1][1] = 0.5
    else:
      let s = sqrt(1.0/3.0)
      im[0][0] = -0.5*s; im[1][1] = -0.5*s; im[2][2] = s
    result[k] = cm(re,im)
proc linear(ms:openArray[C],cs:openArray[float64]):C =
  for k in 0..<ms.len: result = result+cs[k]*ms[k]
proc bridge(m:C, ps:array[8,array[8,C]]):R[8] =
  for i in 0..<8:
    for j in 0..<8: result[i][j] = -tp(ps[i][j],m)
proc invLog[n:static int](a:R[n]):tuple[inv:R[n],log:float] =
  # Pivoted Gauss-Jordan is independent of the production unpivoted LU.
  var x=a
  var sign=1
  result.inv=eye[n]()
  for k in 0..<n:
    var p=k
    for i in k+1..<n:
      if abs(x[i][k])>abs(x[p][k]): p=i
    if p!=k:
      swap(x[p],x[k]); swap(result.inv[p],result.inv[k]); sign = -sign
    let piv=x[k][k]
    doAssert piv!=0
    if piv<0: sign = -sign
    result.log += ln(abs(piv))
    for j in 0..<n:
      x[k][j] /= piv
      result.inv[k][j] /= piv
    for i in 0..<n:
      if i==k: continue
      let c=x[i][k]
      for j in 0..<n:
        x[i][j] -= c*x[k][j]
        result.inv[i][j] -= c*result.inv[k][j]
  doAssert sign>0

proc phiRef(x:R[8], dirs:openArray[R[8]], order=0, scale=0):tuple[p:R[8],dp:seq[R[8]]] =
  # Dense powers and forward jets. order=0 uses converged exp/Phi series.
  var sc=scale
  if order==0:
    sc=0
    var size=normRef(x)
    while size>0.25: size *= 0.5; inc sc
  let h=1.0/float(1 shl sc)
  let y=h*x
  var t=eye[8]()
  var e=t
  var dt=newSeq[R[8]](dirs.len)
  var de=newSeq[R[8]](dirs.len)
  result.p=t
  result.dp.newSeq(dirs.len)
  for k in 1..(if order==0: 256 else: order):
    let den=float(if order==0: k else: k+1)
    var err=0.0
    for d in 0..<dirs.len:
      dt[d]=(1.0/den)*(mm(dt[d],y)+mm(t,h*dirs[d]))
      err=max(err,normRef(dt[d]))
    t=(1.0/den)*mm(t,y)
    if order==0:
      e=e+t
      result.p=result.p+(1.0/float(k+1))*t
      for d in 0..<dirs.len:
        de[d]=de[d]+dt[d]
        result.dp[d]=result.dp[d]+(1.0/float(k+1))*dt[d]
      if k>=4 and max(err,normRef(t))<1e-18: break
      doAssert k<256, "analytic series did not converge"
    else:
      result.p=result.p+t
      for d in 0..<dirs.len: result.dp[d]=result.dp[d]+dt[d]
  for j in 0..<sc:
    let p=result.p
    if order==0:
      for d in 0..<dirs.len:
        result.dp[d]=0.5*(result.dp[d]+mm(de[d],p)+mm(e,result.dp[d]))
        de[d]=mm(de[d],e)+mm(e,de[d])
      result.p=0.5*(p+mm(e,p))
      e=mm(e,e)
    else:
      let c=1.0/float(1 shl (sc+1-j))
      let pp=mm(p,p)
      for d in 0..<dirs.len:
        let dp=result.dp[d]
        result.dp[d]=dp+c*(mm(dirs[d],pp)+mm(x,mm(dp,p)+mm(p,dp)))
      result.p=p+c*mm(x,pp)

proc branch(f:C):int =
  var n2=0.5*dotRef(f,f)
  while n2>1.0/16.0: n2 *= 0.25; inc result
proc expRef(f:C, dirs:openArray[C], ns = -1):tuple[e:C,ds:seq[C],ns:int] =
  result.ns=if ns<0: branch(f) else: ns
  let h=1.0/float(1 shl result.ns)
  let y=h*f
  var t=eye[6]()
  var dt=newSeq[C](dirs.len)
  result.e=t
  result.ds.newSeq(dirs.len)
  for k in 1..12:
    for d in 0..<dirs.len: dt[d]=(1.0/float(k))*(mm(dt[d],y)+mm(t,h*dirs[d]))
    t=(1.0/float(k))*mm(t,y)
    result.e=result.e+t
    for d in 0..<dirs.len: result.ds[d]=result.ds[d]+dt[d]
  for _ in 0..<result.ns:
    for d in 0..<dirs.len: result.ds[d]=mm(result.e,result.ds[d])+mm(result.ds[d],result.e)
    result.e=mm(result.e,result.e)
proc expCase(name:string, f,dm,alpha,cot:C, ns = -1):ExpRef =
  let r=expRef(f,[tah(dm),alpha],ns)
  ExpRef(name:name,f:toMat(f),dm:toMat(dm),alpha:toMat(alpha),c:toMat(cot),
    e:toMat(r.e),de:toMat(r.ds[0]),da:toMat(r.ds[1]),branch:r.ns,
    checkBranch:name.startsWith("threshold") or name in ["zero","tiny","generic012"],
    field:0.5*dotRef(cot,r.ds[0]),alphaPair:0.5*dotRef(cot,r.ds[1]))

let gs=generators()
var products,comms:array[8,array[8,C]]
for i in 0..<8:
  for j in 0..<8: products[i][j]=mm(gs[i],gs[j])
for i in 0..<8:
  for j in 0..<8: comms[i][j]=products[i][j]-products[j][i]

proc caseRef(name:string,m,dm,cot:C,b:R[8],v:array[8,float64],orders:openArray[int]):CaseRef =
  let f=tah(m)
  let x=bridge(m,comms)
  let d=bridge(m,products)
  let dx=bridge(dm,comms)
  let dd=bridge(dm,products)
  let pa=phiRef(x,[dx])
  let ka=eye[8]()+mm(pa.p,d)
  let (ika,loga)=invLog(ka)
  var edirs:array[8,C]
  for i in 0..<8: edirs[i]=tah(mm(gs[i],m))
  let ex=expRef(f,edirs)
  var jj:R[8]
  for j in 0..<8:
    let z=mm(ex.ds[j]+mm(ex.e,gs[j]),trans(ex.e))
    for i in 0..<8: jj[i][j] = -tp(gs[i],z)
  let (_,selected)=invLog(jj)
  result=CaseRef(name:name,m:toMat(m),dm:toMat(dm),c:toMat(cot),x:toAdj(x),dx:toAdj(dx),
    d:toAdj(d),b:toAdj(b),v:toVec(v),fNorm:cnorm(f),dNorm:normRef(d),analyticPhi:toAdj(pa.p),
    analyticLog:loga,analyticDLog:tp(ika,mm(pa.dp[0],d)+mm(pa.p,dd)),
    analyticInvNorm:normRef(ika),analyticCond:normRef(ka)*normRef(ika),selectedLog:selected)
  var dirs:array[9,R[8]]
  for i in 0..<8: dirs[i]=bridge(gs[i],comms)
  dirs[8]=dx
  for order in orders:
    for scale in [0,5]:
      let pr=phiRef(x,dirs,order,scale)
      let k=eye[8]()+mm(pr.p,d)
      let (ik,log)=invLog(k)
      var grad:C
      for i in 0..<3:
        for j in 0..<3:
          for imag in [false,true]:
            var re,im:R[3]
            if imag: im[i][j]=1
            else: re[i][j]=1
            let db=cm(re,im)
            var dp:R[8]
            for a in 0..<8: dp=dp+(-tp(gs[a],db))*pr.dp[a]
            let v=tp(ik,mm(dp,d)+mm(pr.p,bridge(db,products)))
            if imag: grad[i+3][j]=v; grad[i][j+3] = -v
            else: grad[i][j]=v; grad[i+3][j+3]=v
      var apply,cv,cp:array[8,float64]
      for i in 0..<8: cv[i] = -tp(gs[i],tah(cot))
      for i in 0..<8:
        for j in 0..<8:
          apply[i] += pr.p[j][i]*v[j]
          cp[i] += pr.p[j][i]*cv[j]
      result.finite.add FiniteRef(order:order,scale:scale,phi:toAdj(pr.p),dphi:toAdj(pr.dp[8]),
        k:toAdj(k),grad:toMat(grad),apply:toVec(apply),pullback:toMat(linear(gs,cp)),log:log,
        dlog:tp(ik,mm(pr.dp[8],d)+mm(pr.p,dd)),pair:dotRef(b,pr.dp[8]),
        invNorm:normRef(ik),cond:normRef(k)*normRef(ik))

proc buildRefs*():References =
  var cs,v:array[8,float64]
  var b:R[8]
  for i in 0..<8:
    cs[i]=float(i+1)/13
    v[i]=float((5*i) mod 9-4)/17
    for j in 0..<8: b[i][j]=float((3*i+5*j) mod 17-8)/19
  var unit=linear(gs,cs)
  unit=(1/cnorm(unit))*unit
  var dr,di,cr,ci:R[3]
  for i in 0..<3:
    for j in 0..<3:
      dr[i][j]=float(1+2*i-j)/17; di[i][j]=float(2-i+3*j)/19
      cr[i][j]=float((2*i+5*j) mod 7-3)/11; ci[i][j]=float((3*i+2*j) mod 9-4)/13
  var dm=cm(dr,di)
  dm=(1/cnorm(dm))*dm
  let cot=cm(cr,ci)
  let h=cm([[0.09,0.013,-0.008],[0.013,0.11,0.017],[-0.008,0.017,0.14]],
           [[0.0,0.021,0.011],[-0.021,0.0,-0.014],[-0.011,0.014,0.0]])
  var zr,rep,ext:R[3]
  rep[0][0]=1; rep[1][1]=1; rep[2][2] = -2
  ext[0][0]=1; ext[1][1] = -1
  let repeated=(1/sqrt(6.0))*cm(zr,rep)
  let extremal=(1/sqrt(2.0))*cm(zr,ext)
  let q=cm(mm([[0.8,-0.6,0.0],[0.6,0.8,0.0],[0.0,0.0,1.0]],
              [[1.0,0.0,0.0],[0.0,12.0/13,-5.0/13],[0.0,5.0/13,12.0/13]]),zr)
  let rotated=mm(mm(q,extremal),trans(q))
  let specs=[("zero",unit,0.0),("tiny",unit,1e-12),("generic012",unit,0.12),
    ("generic1",unit,1.0),("generic2",unit,2.0),("generic4",unit,4.0),("generic8",unit,8.0),
    ("repeated4",repeated,4.0),("repeated8",repeated,8.0),("extremal8",extremal,8.0),
    ("rotated4",rotated,4.0),("rotated8",rotated,8.0)]
  for (name,u,amp) in specs:
    let m=amp*u+h
    let orders=if name=="generic012": @[1,3,7,11,13] else: @[13]
    result.cases.add caseRef(name,m,dm,cot,b,v,orders)
    result.exp.add expCase(name,tah(m),dm,u,cot)
  var thresholds:seq[tuple[name:string,f:C]]
  for amp in [0.25,0.5,1.0,2.0,4.0]:
    for side in [-1,1]:
      let f=(amp*(1+float(side)*1e-8))*unit
      let name="threshold" & $amp & (if side<0: "below" else: "above")
      thresholds.add (name,f)
      result.exp.add expCase(name,f,dm,unit,cot)
  var ns=0
  for i in 0..<4: ns=max(ns,branch(thresholds[i].f))
  for i in 0..<4: result.simd.add expCase(thresholds[i].name,thresholds[i].f,dm,unit,cot,ns)

