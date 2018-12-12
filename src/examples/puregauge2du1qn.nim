import qex
import gauge, physics/qcdTypes
import mdevolve
import os, strutils, algorithm, macros
import random  # for seeding our RNGs

const CHECK = true
const CHECKLBFGS = false
const CHECKLBFGSMOM = true
const CHECKLBFGSEIGEN = true

when CHECKLBFGSEIGEN:
  when not defined primmeDir:
    static: error "Must declare primmeDir to CHECKLBFGSEIGEN"
  when not defined lapackLib:
    static: error "Must declare lapackLib to CHECKLBFGSEIGEN"
  import eigens/qexPrimmeInternal
  import primme

template checkgc =
  let fl = instantiationInfo(fullPaths=true)
  echo "CHECK GC: ",fl.filename,":",fl.line
  GC_fullCollect()
  echo GC_getStatistics()
  #dumpNumberOfInstances()  # remember to -d:nimTypeNames

proc norm2(y:seq[DLatticeColorMatrixV]):auto =
  var r {.noinit.}: type(y[0].norm2)
  r := 0.0
  for mu in 0..<y.len:
    r += y[mu].norm2
  r

proc dot(x,y:seq[DLatticeColorMatrixV]):auto =
  var r {.noinit.}: type(x[0].dot y[0])
  r := 0.0
  for mu in 0..<x.len:
    r += x[mu].dot y[mu]
  r

proc redot(x,y:seq[DLatticeColorMatrixV]):auto =
  var r {.noinit.}: type(x[0].redot y[0])
  r := 0.0
  for mu in 0..<x.len:
    r += x[mu].redot y[mu]
  r

proc setup(defaultLat:openarray[int],
    nstream:int = 9,
    seed:int64 = int64(11^11)):auto =
  var lat:seq[int]
  let pc = paramCount()
  if pc > 0 and paramStr(1).isDigit:
    lat = @[]
    for i in 1..pc:
      if not paramStr(i).isDigit: break
      lat.add paramStr(i).parseInt
  else:
    lat = @defaultLat
  let
    lo = lat.newLayout
  var
    seedrng = initRand seed
    gs = newseq[type(lo.newGauge)](nstream)
  let
    fseed = seedrng.next
    gseed = seedrng.next
  echo "fieldSeed = ",fseed
  echo "globalSeed = ",gseed
  for i in 0..<nstream: gs[i] = lo.newGauge
  var
    r = newRNGField(RngMilc6, lo, fseed)
    R:RngMilc6
  R.seed(gseed, 0)
  return (gs, r, R)

proc topo2DU1(g:array or seq):float =
  tic()
  const nc = g[0][0].nrows
  let
    lo = g[0].l
    nd = lo.nDim
    t = newTransporters(g, g[0], 1)
  var p = 0.0
  toc("topo2DU1 setup")
  threads:
    tic()
    var tp:type(atan2(g[0][0][0,0].im, g[0][0][0,0].re))
    for mu in 1..<nd:
      for nu in 0..<mu:
        let tpl = (t[mu]^*g[nu]) * (t[nu]^*g[mu]).adj
        for i in tpl:
          tp += atan2(tpl[i][0,0].im, tpl[i][0,0].re)
    var v = tp.simdSum
    v.threadRankSum
    threadSingle: p += v
    toc("topo2DU1 work")
  toc("topo2DU1 threads")
  p/TAU

proc maxTreeFix(f:seq, val:float, fixextra:bool) =
  ## Set link to `val`, for those links on the maximal tree.
  let
    nd = f.len
    lo = f[0].l
    lat = lo.physGeom
  threads:
    var co = newseq[cint](nd)
    for i in 0..<nd:
      for j in lo.sites:
        lo.coord(co,(lo.myRank,j))
        var zeroafter = true
        for k in i+1 .. nd-1: zeroafter = zeroafter and co[k] == 0
        if zeroafter:
          if co[i] < lat[i]-1:
            f[i]{j} := val
          elif fixextra:
            var zerobefore = true
            for k in 0..<i: zerobefore = zerobefore and co[k] == 0
            if zerobefore: f[i]{j} := val

type LBFGS[F] = object
  ys,gamma:seq[float]
  y,s,v:seq[F]  ## v_k = y_k + beta_k G_k-1 s_k
  sortedix:seq[int]
  p:int  ## Index for a new addition, or 1 plus the end of the ring.  Temporary storage at p-1.
  sortedlen:int
  kappa:float  ## minimum eigen value of H to regularize the zero modes.
  invH0:float
  first:bool

proc initLBFGS(lo:Layout, n:int, h0 = 1.0, kappa:float = 0.0):auto =
  type F = type(lo.newGauge)
  var r {.noinit.} :LBFGS[F]
  r.p = 0
  r.first = true
  r.ys.newseq(n)
  r.y.newseq(n)
  r.s.newseq(n)
  r.sortedlen = 0
  r.sortedix.newseq(n-2)  # N states -> N differences ->  N-2 excluding 1 state
  r.v.newseq(n-2)  # ditto
  r.gamma.newseq(n-2)
  r.invH0 = 1.0 / h0
  r.kappa = kappa
  for i in 0..<n:
    r.y[i] = lo.newGauge
    r.s[i] = lo.newGauge
  for i in 0..<n-2:
    r.v[i] = lo.newGauge
  r

proc add[F](o:var LBFGS[F]; x,f:F) =
  var n = o.p
  let
    k = o.kappa
    k1 = 1.0 - k
  if o.first:
    o.first = false
  else:
    # save new s = ln(x_old x.adj), y = k s + (1-k) (f_old-f)
    let n1 = if n == 0: o.ys.len - 1 else: n - 1
    o.ys[n1] = 0.0
    for mu in 0..<x.len:
      for i in x[mu]:
        when CHECKLBFGS:
          var t = o.s[n1][mu][i] - x[mu][i]
        else:
          var t = ln(o.s[n1][mu][i] * x[mu][i].adj)
          t.projectTAH t
          t := t.im
        when CHECK and false:
          if i == 0:
            echo "For mu = ",mu
            let
              si = o.s[n1][mu][i]
              xi = x[mu][i]
              xia = xi.adj
            echo "xi_old: ",si
            echo "xi: ",xi
            echo "xi.adj: ",xia
            echo "ln(xi_old*xi.adj): ",ln(si * xia)
            echo "==t: ",t
        o.s[n1][mu][i] := t
        when CHECKLBFGS:
          var r = o.y[n1][mu][i] - f[mu][i]
        else:
          var r = (o.y[n1][mu][i] - f[mu][i]).im
        r = k * t + k1 * r  # FIXME: this is useless
        when CHECK and false:
          if i == 0:
            let
              yi = o.y[n1][mu][i]
              fi = f[mu][i]
            echo "fi_old: ",yi
            echo "fi: ",fi
            echo "fi-fi_old: ",fi-yi
            echo "==r: ",r
        o.y[n1][mu][i] := r
      o.ys[n1] += o.y[n1][mu].redot o.s[n1][mu]
    when CHECK and false:
      let ysd = o.y[n1].dot o.s[n1]
      echo "ysd: ",ysd
      echo "ys[",n1,"]: ",o.ys[n1]
  for mu in 0..<x.len:
    o.s[n][mu] := x[mu]  # TODO: we can elide this extra copy
    o.y[n][mu] := f[mu]
  inc n
  if n == o.ys.len: n = 0
  o.p = n

proc reverseadd[F](o:var LBFGS[F]; x,f:F) =
  # Eg. 4 states, after one forward pass, we have saved
  # 0'    1-2  2-3  3-0',  p = 1
  # 0'-1'  1'    2-3  3-0',  p = 2
  # 0'-1'  1'-2'  2'    3-0',  p = 3
  # Now revert the order, and we need
  # 2'    1'-2'  0'-1'  0',  with p = 0
  # Then after the next (3'->3), we will have subsequently
  # 3    1'-2'  0'-1'  0'-3, with p = 1
  # 3-2  2    0'-1'  0'-3,  with p = 2
  # 3-2  2-1  1    0'-3, with p = 3
  # Note that the commutation of the difference does not matter,
  # as long as y and s are consistent.
  let p = o.p
  if p != o.ys.len-1:
    echo "reverseadd: p is not at the end."
    qexAbort()
  for i in 0..<(p div 2):
    let t = o.ys[i]
    o.ys[i] = o.ys[p-1-i]
    o.ys[p-1-i] = t
    for mu in 0..<o.y[p].len:
      for e in o.y[i][mu]:
        let t = o.y[i][mu][e]
        o.y[i][mu][e] := o.y[p-1-i][mu][e]
        o.y[p-1-i][mu][e] := t
        let r = o.s[i][mu][e]
        o.s[i][mu][e] := o.s[p-1-i][mu][e]
        o.s[p-1-i][mu][e] := r
  for mu in 0..<x.len:
    o.s[p][mu] := x[mu]
    o.y[p][mu] := f[mu]
  # o.ys[p] doesn't matter.
  o.p = 0

proc A[F](o:LBFGS[F], k:int, z:var F) =
  ## A_k z  ->  z  where  H_k = A_k A_k^\dag
  # A0 z -> z
  let a0 = 1.0 / sqrt o.invH0
  for mu in 0..<z.len: z[mu] *= a0
  for j in 0..k:
    let i = o.sortedix[j]
    let sz = (-o.gamma[j]) * o.s[i].dot z
    for mu in 0..<z.len:
      for e in z[mu]:
        let t = z[mu][e] + sz * o.v[j][mu][e]
        z[mu][e] := t
proc Adag[F](o:LBFGS[F], k:int, z:var F) =
  ## A_k^\dag z  ->  z  where  H_k = A_k A_k^\dag
  for j in countdown(k,0):
    let i = o.sortedix[j]
    let vz = (-o.gamma[j]) * o.v[j].dot z
    for mu in 0..<z.len:
      for e in z[mu]:
        let t = z[mu][e] + vz * o.s[i][mu][e]
        z[mu][e] := t
  # A0d z -> z
  let a0 = 1.0 / sqrt o.invH0
  for mu in 0..<z.len: z[mu] *= a0
proc H[F](o:LBFGS[F], k:int, r:var F) =
  when CHECKLBFGS and false: echo "r: ",r
  o.Adag(k, r)
  when CHECKLBFGS and false: echo "Adr: ",r
  o.A(k, r)
  when CHECKLBFGS and false: echo "AAdr: ",r
proc H[F](r:var F, o:LBFGS[F], k:int, z:F) =
  for mu in 0..<z.len: r[mu] := z[mu]
  o.H(k,r)
proc H[F](o:LBFGS[F], r:var F) =
  let n = o.sortedlen
  o.H(n-1, r)
proc H[F](r:var F, o:LBFGS[F], z:F) =
  let n = o.sortedlen
  r.H(o, n-1, z)
proc sqrtH[F](o:LBFGS[F]; x:var F) =
  ## A x  ->  x  where  H = A A^\dag
  let n = o.sortedlen
  o.A(n-1, x)

proc B[F](o:LBFGS[F], k:int, z:var F) =
  ## B_k z  ->  z  where  H_k^-1 = B_k B_k^\dag
  # B0 z -> z
  let b0 = sqrt o.invH0
  for mu in 0..<z.len: z[mu] *= b0
  for j in 0..k:
    let i = o.sortedix[j]
    let vz = (o.v[j].dot z) / (-o.ys[i])
    for mu in 0..<z.len:
      for e in z[mu]:
        let t = z[mu][e] + vz * o.s[i][mu][e]
        z[mu][e] := t
proc Bdag[F](o:LBFGS[F], k:int, z:var F) =
  ## A_k^\dag z  ->  z  where  H_k^-1 = B_k B_k^\dag
  for j in countdown(k,0):
    let i = o.sortedix[j]
    let sz = (o.s[i].dot z) / (-o.ys[i])
    for mu in 0..<z.len:
      for e in z[mu]:
        let t = z[mu][e] + sz * o.v[j][mu][e]
        z[mu][e] := t
  # B0d z -> z
  let b0 = sqrt o.invH0
  for mu in 0..<z.len: z[mu] *= b0
proc invH[F](o:LBFGS[F], k:int, r:var F) =
  when CHECKLBFGS and false: echo "r: ",r
  o.Bdag(k, r)
  when CHECKLBFGS and false: echo "Bdr: ",r
  o.B(k, r)
  when CHECKLBFGS and false: echo "BBdr: ",r
proc invH[F](r:var F, o:LBFGS[F], k:int, z:F) =
  for mu in 0..<z.len: r[mu] := z[mu]
  o.invH(k,r)
proc invH[F](o:LBFGS[F], r:var F) =
  let n = o.sortedlen
  o.invH(n-1, r)
proc invH[F](r:var F, o:LBFGS[F], z:F) =
  let n = o.sortedlen
  r.invH(o, n-1, z)
proc sqrtInvH[F](o:LBFGS[F]; x:var F) =
  ## B x  ->  x  where  H^-1 = B B^\dag
  let n = o.sortedlen
  o.B(n-1, x)

proc invHl[F](o:LBFGS[F]; x:var F) =
  ## The basic lBGFS 2-loop algorithm.
  ## H^-1 y  ->  x
  let n = o.sortedlen
  type C = type(x.dot x)
  var alpha = newseq[C](n)
  for j in countdown(n-1,0):
    let k = o.sortedix[j]
    let a = (o.s[k].dot x) / o.ys[k]
    alpha[j] := a
    when CHECK and false:
      echo "invH ",j," (",k,")"
      echo "x.norm2: ",x.norm2
      echo "alpha: ",a
      echo "a^2*o.y[",k,"].norm2: ",(a*a*o.y[k].norm2)
    for mu in 0..<x.len:
      for e in x[mu]:
        let ay = a * o.y[k][mu][e]
        let t = x[mu][e] - ay
        when CHECK and false:
          if e == 0:
            echo "x[",mu,"][0]: ",x[mu][e]
            echo "a*o.y[",k,"][",mu,"][0]: ",ay
            echo "t: ",t
        x[mu][e] := t
    when CHECK and false:
      echo "(1-rys)x.norm2: ",x.norm2
  # Initial inverse Hessian:  x *= YS/YY
  for mu in 0..<x.len: x[mu] *= o.invH0
  for j in 0..<n:
    let k = o.sortedix[j]
    let r = (o.y[k].dot x) / o.ys[k]
    let beta = alpha[j] - r
    when CHECK and false:
      echo "invH2 ",j," (",k,")"
      echo "alpha: ",alpha[j]
      echo "o.y[",k,"].dot(x)/o.ys[",k,"] ",r
      echo "beta: ",beta
      echo "x.norm2: ",x.norm2
    for mu in 0..<x.len:
      for e in x[mu]:
        let bs = beta * o.s[k][mu][e]
        let t = x[mu][e] + bs
        when CHECK and false:
          if e == 0:
            echo "x[",mu,"][0]: ",x[mu][e]
            echo "beta*o.s[",k,"][",mu,"][0]: ",bs
            echo "t: ",t
        x[mu][e] := t
    when CHECK and false:
      echo "(1+bs)x.norm2: ",x.norm2
proc invHl[F](x:var F; o:LBFGS[F]; y:F) =
  ## H^-1 y  ->  x
  for mu in 0..<x.len: x[mu] := y[mu]
  o.invHl x

when CHECKLBFGSEIGEN:
  template genMatvec(op:untyped):untyped =
    proc `lbfgs op Matvec`[MI](x:pointer, ldx:ptr PRIMME_INT, y:pointer, ldy:ptr PRIMME_INT,
        blocksize:ptr cint, primme:ptr primme_params, err:ptr cint) {.noconv.} =
      var
        x = asarray[PRIMME_COMPLEX_DOUBLE] x
        dx = ldx[]
        y = asarray[PRIMME_COMPLEX_DOUBLE] y
        dy = ldy[]
      for i in 0..<blocksize[]:
        let
          xp = x[i*dx].addr         # Input vector
          yp = y[i*dy].addr         # Output
          m = cast[ptr MI](primme.matrix)
        var g = m.tmpgauge
        g.fromPrimmeArrayGauge xp
        m.lbfgs.op g
        g.toPrimmeArrayGauge yp
      err[] = 0
  genMatvec(invH)
  genMatvec(invHl)
  genMatvec(H)
  genMatvec(sqrtH)

proc prep[F](o:var LBFGS[F], cutoff = 0.0, reduce = 0) =
  ## For differences y and s, the item k has the difference between k and k+1.
  ## Before using the approximation, o.p points to the current updating stream number.
  ## We sort according to ys, excluding item o.p and o.p-1, so we don't depend on ourselves.
  let e0 = o.p
  var e1 = e0 - 1
  if e1 < 0: e1 = o.ys.len - 1
  o.sortedlen = o.sortedix.len
  var ix = 0
  for i in 0..<o.sortedix.len:
    if ix == e1: ix += 2
    elif ix == e0: ix += 1
    o.sortedix[i] = ix
    inc ix
  let ys = o.ys
  proc cmpys(i,j:int):int =
    let a = ys[i]
    let b = ys[j]
    result = cmp(a,b)
    #if a <= 0 or b <= 0:
    #  result = -result
  o.sortedix.sort(cmpys, SortOrder.Descending)
  var unset = true
  for i in 0..<o.sortedix.len:
    let k = o.sortedix[i]
    let yy = o.y[k].norm2
    let ss = o.s[k].norm2
    echo i," ix ",k," ys ",o.ys[k]," ys/yy ",(o.ys[k]/yy)," yy ",yy," ss ",ss
    if unset and o.ys[k] <= cutoff:  # FIXME how to tune this?
      o.sortedlen = i
      unset = false
  if reduce > 0:
    echo "lBGFS nmem reduce: ",reduce
    o.sortedlen -= reduce
    if o.sortedlen < 0: o.sortedlen = 0
  echo "lBFGS nmem = ",o.sortedlen
  #if o.sortedlen > 0 and false:
  #  #let k = o.sortedix[o.sortedlen-1]
  #  let k = o.sortedix[0]  # FIXME: Does not help stability
  #  let yy = o.y[k].norm2
  #  when CHECK and false:
  #    #let yr = o.y[k].redot o.y[k]
  #    #let yd = o.y[k].dot o.y[k]
  #    #echo "y.dot y: ",yd
  #    #echo "y.redot y: ",yr
  #    echo "y.norm2: ",yy
  #  o.invH0 = o.ys[k] / yy
  #else:
  #  o.invH0 = 1.0
  var gs = o.s[0].newOneOf
  for j in 0..<o.sortedlen:  # Compute v_k and gamma_k
    # FIXME maybe save more vectors?  possibly linear time?
    let k = o.sortedix[j]
    gs.H(o, j-1, o.s[k])  # G_k-1 s_k  (j is the actual order number)
    let alpha = 1.0 / o.s[k].redot gs
    when CHECK and false:
      let alphad = 1.0 / o.s[k].dot gs
      echo "alpha: ",alpha
      echo "alphad: ",alphad
    let beta = -sqrt(alpha * o.ys[k])  # FIXME: negative sign more stable?
    o.gamma[j] = beta / o.ys[k]
    when CHECK:
      echo j," alpha: ",alpha," beta: ",beta," gamma: ",o.gamma[j]
    for mu in 0..<gs.len:
      for e in gs[mu]:
        let t = o.y[k][mu][e] + beta * gs[mu][e]
        o.v[j][mu][e] := t
  when CHECKLBFGSEIGEN:
    type MatrixInfo = object
      lbfgs: type(o)
      tmpgauge: type(o.y[0])
    const nc = o.y[0][0][0].nrows
    let
      lo = o.y[0][0].l
      nd = lo.nDim
    proc lbfgseigen(o:LBFGS[F],t:primme_target) =
      var p = primme_initialize()
      p.n = nc*nc*nd*lo.physVol
      p.matrixMatvec = lbfgsHMatvec[MatrixInfo]
      p.numEvals = min(o.y.len, p.n div 2).cint
      p.target = t
      p.eps = 1e-13
      p.numProcs = nRanks.cint
      p.procId = myRank.cint
      p.nLocal = nc*nc*nd*lo.nSites
      p.globalSumReal = sumReal[primme_params]
      p.printLevel = 1.cint  # 3.cint
      block primmeSetMethod:
        let ret = p.set_method PRIMME_DEFAULT_METHOD
        if 0 != ret:
          echo "ERROR: set_method returned with nonzero exit status: ", ret
          qexAbort()
      block prepare:
        let ret = zprimme(nil,nil,nil,p.addr)
        if 1 != ret:
          echo "Error: zprimme(nil) returned with exit status: ", ret
          qexAbort()
      var
        vecs = newAlignedMemU[complex[float]]int(p.numEvals*p.nLocal)
        vals = newseq[float]p.numEvals
        rnorms = newseq[float]p.numEvals
        intWork = newAlignedMemU[char]p.intWorkSize
        realWork = newAlignedMemU[char]p.realWorkSize
      var mi = MatrixInfo(lbfgs: o)
      mi.tmpgauge = newOneOf o.y[0]
      p.intWork = cast[ptr cint](intWork.data)
      p.realWork = realWork.data
      p.matrix = mi.addr
      if myRank == 0: p.display_params
      block run:
        let ret = p.run(vals, asarray[PRIMME_COMPLEX_DOUBLE](vecs.data)[], rnorms)
        if ret != 0:
          echo "Error: primme returned with nonzero exit status: ", ret
          qexAbort()
      echo "Neigens    : ",p.initSize
      echo "Iterations : ",p.stats.numOuterIterations
      echo "Restarts   : ",p.stats.numRestarts
      echo "Matvecs    : ",p.stats.numMatvecs
      echo "Preconds   : ",p.stats.numPreconds
      echo "GlobalSums : ",p.stats.numGlobalSum
      echo "VGlobalSum : ",p.stats.volumeGlobalSum
      echo "OrthoIProd : ",p.stats.numOrthoInnerProds
      echo "ElapsedT   : ",p.stats.elapsedTime
      echo "MatvecT    : ",p.stats.timeMatvec
      echo "PrecondT   : ",p.stats.timePrecond
      echo "OrthoT     : ",p.stats.timeOrtho
      echo "GlobalSumT : ",p.stats.timeGlobalSum
      echo "EstMinEv   : ",p.stats.estimateMinEVal
      echo "EstMaxEv   : ",p.stats.estimateMaxEVal
      echo "EstMaxSv   : ",p.stats.estimateLargestSVal
      echo "EstResid   : ",p.stats.estimateResidualError
      echo "MaxConvTol : ",p.stats.maxConvTol
      if p.locking != 0 and p.intWork != nil and p.intWork[] == 1:
        echo "\nA locking problem has occurred."
        echo "Some eigenpairs do not have a residual norm less than the tolerance."
        echo "However, the subspace of evecs is accurate to the required tolerance."
      case p.dynamicMethodSwitch:
      of -1: echo "Recommended method for next run: DEFAULT_MIN_MATVECS"
      of -2: echo "Recommended method for next run: DEFAULT_MIN_TIME"
      of -3: echo "Recommended method for next run: DYNAMIC (close call)"
      else: discard
      for i in 0..<p.numEvals:
        echo i," ev ",vals[i]," rnorm ",rnorms[i]
    o.lbfgseigen(primme_smallest)
    o.lbfgseigen(primme_largest)
    #qexExit 0

when CHECKLBFGS:
    const cklbfgsnd = 64
    let nseq = 65536
    type V = array[cklbfgsnd,float]
    type M = array[cklbfgsnd,V]
    type D = distinct V
    template len(x:D):int = V(x).len
    template `[]`(x:D,i:untyped):untyped = V(x)[i]
    template `[]=`(x:D,i:untyped,y:untyped):untyped = V(x)[i] = y
    template `$`(x:D):string = $V(x)
    iterator items(x:D):int =
      for i in 0..<x.len: yield i
    proc `*`(c:float,x:D):D =
      for i in x:
        result[i] = c * x[i]
    proc `+`(x,y:D):D =
      for i in x:
        result[i] = x[i] + y[i]
    proc `-`(x,y:D):D =
      for i in x:
        result[i] = x[i] - y[i]
    proc dot(x,y:D):float =
      for i in x:
        result += x[i]*y[i]
    proc `*=`(x:var D,y:float) =
      for i in x:
        x[i] *= y
    proc `-=`(x:var D,y:D) =
      for i in x:
        x[i] -= y[i]
    template `:=`(x,y:float) = x = y
    template `:=`(x,y:D) =
      for i in x:
        x[i] = y[i]
    template redot(x,y:D):float = x.dot y
    template norm2(x:D):float = x.dot x
    proc newOneOf(x:D):D = x
    type F = seq[D]
    proc dot(x,y:F):float =
      for i in 0..<x.len:
        result += x[i].dot y[i]
    template redot(x,y:F):float = x.dot y
    template norm2(y:F):float =
      let x = y
      x.dot x
    proc `*`(m:M,y:D):D =
      for j in 0..<m.len:
        for i in 0..<m.len:
          result[i] += m[i][j] * y[j]
    var cklbfgsr {.noinit.} :LBFGS[F]
    cklbfgsr.p = 0
    cklbfgsr.first = true
    cklbfgsr.ys.newseq(nseq)
    cklbfgsr.y.newseq(nseq)
    cklbfgsr.s.newseq(nseq)
    cklbfgsr.sortedlen = 0
    cklbfgsr.sortedix.newseq(nseq-2)  # N states -> N differences ->  N-2 excluding 1 state
    cklbfgsr.v.newseq(nseq-2)  # ditto
    cklbfgsr.gamma.newseq(nseq-2)
    cklbfgsr.invH0 = 1.0
    for i in 0..<nseq:
      cklbfgsr.y[i].newseq(1)
      cklbfgsr.s[i].newseq(1)
    for i in 0..<nseq-2:
      cklbfgsr.v[i].newseq(1)
    var rng = initRand 5^5
    var m:M
    for i in 0..<m.len-1:
      m[i][i] = float((i+1)*(i+1))
      m[i][i+1] = 1
      m[i+1][i] = 1
    m[m.len-1][m.len-1] = 1
    for i in 0..<nseq:
      var x:D
      for k in 0..<cklbfgsnd: x[k] = rng.rand 1.0
      #echo "x: ",x
      let y = m * x
      #echo "y: ",y
      cklbfgsr.add(x = @[x], f = @[y])
    cklbfgsr.prep
    var x:D
    for i in x: x[i] = 1
    echo "test x: ",x
    let y = m * x
    echo "test y: ",y
    var hx = @[x]
    cklbfgsr.H hx
    echo "hx: ",hx
    echo "hx err: ",norm2(hx[0]-y)
    var hinvhx = hx
    cklbfgsr.invH hinvhx
    echo "hinvhx: ",hinvhx
    echo "hinvhx err: ",norm2(hinvhx[0]-x)
    var hlinvhx = hx
    cklbfgsr.invHl hlinvhx
    echo "hlinvhx: ",hlinvhx
    echo "hlinvhx err: ",norm2(hlinvhx[0]-x)
    var hinvy = @[y]
    cklbfgsr.invH hinvy
    echo "hinvy: ",hinvy
    echo "hinvy err: ",norm2(hinvy[0]-x)
    var hlinvy = @[y]
    cklbfgsr.invHl hlinvy
    echo "hlinvy: ",hlinvy
    echo "hlinvy err: ",norm2(hlinvy[0]-x)
    var hsqrtx = @[x]
    cklbfgsr.sqrtH hsqrtx
    echo "hsqrtx: ",hsqrtx
    qexExit 0

qexinit()
threads: echo "thread ",threadNum," / ",numThreads

let
  beta = floatParam("beta", 5.0)
  trajs = intParam("trajs", 512)
  qnbegin = intParam("qnbegin",64)
  tau = floatParam("tau", 2.0)
  steps = intParam("steps", 10)
  qntau = floatParam("qntau", tau * sqrt(2*beta))  # sqrt diagonal term of the Hessian
  qnsteps = intParam("qnsteps", steps)
  qnyscut = floatParam("qnyscut",9)
  seed = intParam("seed", 11^11)
  gfix = intParam("gfix", 0).bool
  gfixextra = intParam("gfixextra", 0).bool
  gfixunit = intParam("gfixunit", 1).bool
  nstream = intParam("nstream", 10)
  kappa = floatParam("kappa", 0)

echo "beta = ",beta
echo "trajs = ",trajs
echo "tau = ",tau
echo "steps = ",steps
echo "qnbegin = ",qnbegin
echo "qntau = ",qntau
echo "qnsteps = ",qnsteps
echo "qnyscut = ",qnyscut
echo "gfix = ",gfix.int
echo "gfixextra = ",gfixextra.int
echo "gfixunit = ",gfixunit.int
echo "nstream = ",nstream
echo "kappa = ",kappa
echo "seed = ",seed

var (gs,r,R) = setup([64,64],nstream,seed)

let
  lo = gs[0][0].l
  lat = lo.physGeom
  nd = lat.len
  gc = GaugeActionCoeffs(plaq:beta)

echo "latsize = ",lat
echo "volume = ",lo.physVol

template getforce(f,g:untyped) =
  f.gaugeforce2(g, gc)
  if gfix: f.maxTreeFix(0.0, gfixextra)

var
  lbfgs = lo.initLBFGS(nstream, h0 = 2*beta, kappa = kappa)
  p = lo.newgauge
  f = lo.newgauge
  g0 = lo.newgauge
  nsnow = 0  # The serial number of the current running HMC.
  forward = true
  hinvp = lo.newgauge

for i in 0..<nstream:
  gs[i].random r
  if gfix and gfixunit: gs[i].maxTreeFix(1.0, gfixextra)
  echo "Initial plaq@gs[",i,"]: ",gs[i].plaq3

proc mdt(t:float) =
  if lbfgs.sortedlen == 0:
    threads:
      for i in 0..<gs[nsnow].len:
        for e in gs[nsnow][i]:
          let etpg = exp((-t)*p[i][e])*gs[nsnow][i][e]
          gs[nsnow][i][e] := etpg
  else:
    hinvp.invH(lbfgs, p)
    when CHECK and false:
      var pn2 = hinvp.norm2
      echo "p.norm2: ",p.norm2,"  hinvp.norm2: ",pn2,"  p.hinv.p.norm2: ",(p.redot hinvp)
      #echo "p[0][0]: ",p[0][0]
      #echo "p[1][0]: ",p[1][0]
      echo "hinvp[0][0]: ",hinvp[0][0]
      #echo "hinvp[1][0]: ",hinvp[1][0]
    #hinvp.projectTAH
    #when CHECK:
    #  let ppn2 = hinvp.norm2
    #  let e = abs(pn2-ppn2)/max(pn2.abs,ppn2.abs)
    #  echo "hinvp.p[0][0]: ",hinvp[0][0]
    #  echo "hinvp.p[1][0]: ",hinvp[1][0]
    #  if e > 1e-12: echo "ProjectTAH changed p.norm2: ",e
    threads:
      for i in 0..<gs[nsnow].len:
        for e in gs[nsnow][i]:
          let etp = exp((-t)*hinvp[i][e])
          when CHECK and false:
            if e == 0:
              echo "etp[",i,"][0]: ",etp
          let etpg = etp*gs[nsnow][i][e]
          gs[nsnow][i][e] := etpg
proc mdv(t:float) =
  f.getforce gs[nsnow]
  threads:
    for i in 0..<f.len:
      for e in f[i]:
        let tf = t*f[i][e]
        p[i][e] += tf


var
  # md = mkLeapfrog(steps = steps, V = mdv, T = mdt)
  # md = mkSW92(steps = steps, V = mdv, T = mdt)
  md = mkOmelyan2MN(steps = steps, V = mdv, T = mdt)
  # md = mkOmelyan4MN4FP(steps = steps, V = mdv, T = mdt)
  # md = mkOmelyan4MN5FV(steps = steps, V = mdv, T = mdt)

for n in 1..trajs:
  if n == qnbegin:
    echo "STARTING QN update"
    md.steps = qnsteps
    for i in 0..<nstream:
      f.getforce gs[i]
      lbfgs.add(x = gs[i], f = f)
    when CHECKLBFGSMOM:
      # Refreshing the momentum and check lbfgs stability.
      var pcheck = newOneOf p
      proc checklbfgsmom =
        var p2 = 0.0
        threads: p.randomTAH r
        if gfix: p.maxTreeFix(0.0, gfixextra)  # FIXME maybe only one is needed
        threads:
          var p2t = 0.0
          for i in 0..<p.len:
            p2t += p[i].norm2
          threadMaster: p2 = p2t
        lbfgs.sqrtH p
        if gfix: p.maxTreeFix(0.0, gfixextra)  # FIXME
        hinvp.invH(lbfgs, p)
        pcheck.invHl(lbfgs, p)
        #threads:
        block:
          var p2t,p2tc,p2nt,hp2 = 0.0
          for i in 0..<p.len:
            p2t += p[i].redot hinvp[i]
            p2tc += p[i].redot pcheck[i]
            p2nt += hinvp[i].norm2
            hp2 += p[i].norm2
          #threadMaster:
          block:
            let e = abs(p2-p2t)/max(p2.abs,p2t.abs)
            let el = abs(p2-p2tc)/max(p2.abs,p2tc.abs)
            echo "SqrtHP.norm2: ",hp2," HinvP.norm2: ",p2nt
            echo "p2: ",p2
            echo "p2_Hinv: ",p2t," err: ",e
            echo "p2_Hinvl: ",p2tc," err: ",el
      for c in countdown(10,0):
        let cut = float(c*lat[0])
        echo "LBFGS cutoff: ",cut
        lbfgs.prep(cutoff = cut)
        for i in 0..<1024: checklbfgsmom()
      echoTimers()
      checkgc()
      qexExit 0

  for ns in 0..<nstream:
    nsnow = if forward: ns else: nstream-1-ns
    echo "Begin traj: ",n," nsNow: ",nsnow," ",(if forward: "forward" else: "backward")
    if n >= qnbegin:
      lbfgs.prep(cutoff = qnyscut * lat[0].float)
      when CHECK and false:
        # Goodness of the Hessian approximation
        let dev = 0.001
        let s0 = gs[nsnow].gaugeAction2 gc
        f.getforce gs[nsnow]
        for i in 0..<g0.len:
          for j in g0[i]:
            let t = dev * ln gs[nsnow][i][j]
            g0[i][j] := t
        let s01 = f.redot g0
        lbfgs.sqrtH g0
        if gfix: g0.maxTreeFix(0.0, gfixextra)
        let s02 = 0.5 * g0.norm2
        let s1a = s0 + s01 + s02
        for i in 0..<g0.len:
          for j in g0[i]:
            let t = exp((1+dev) * ln gs[nsnow][i][j])
            g0[i][j] := t
        if gfix: g0.maxTreeFix(1.0, gfixextra)
        let s1 = g0.gaugeAction2 gc
        echo "s0: ",s0," s01: ",s01," s02: ",s02
        echo "s1a: ",s1a
        echo "s1: ",s1
        echo "error: ",(s1-s1a)
        echo "s1-(s0+s01): ",(s1-(s0+s01))

    var p2 = 0.0
    threads: p.randomTAH r
    if gfix: p.maxTreeFix(0.0, gfixextra)  # FIXME maybe only one is needed
    if CHECK or n < qnbegin:
      threads:
        var p2t = 0.0
        for i in 0..<p.len:
          p2t += p[i].norm2
        threadMaster: p2 = p2t

    if n >= qnbegin:
      lbfgs.sqrtH p
      if gfix: p.maxTreeFix(0.0, gfixextra)  # FIXME
      hinvp.invH(lbfgs, p)
      threads:
        var p2t = 0.0
        when CHECK:
          var p2nt,hp2 = 0.0
        for i in 0..<p.len:
          p2t += p[i].redot hinvp[i]
          when CHECK:
            p2nt += hinvp[i].norm2
            hp2 += p[i].norm2
        threadMaster:
          when CHECK:
            echo "SqrtHP.norm2: ",hp2
            echo "HinvP.norm2: ",p2nt
            echo "p2 simple: ",p2
            echo "p2 with_Hinv: ",p2t
            var e = abs(p2-p2t)/max(p2.abs,p2t.abs)
            if e > 1e-12:
              echo "ERROR: Failed lbfgs consistency check: e = ",e
          p2 = p2t

    threads:
      for i in 0..<g0.len:
        g0[i] := gs[nsnow][i]

    let
      ga0 = g0.gaugeAction2 gc
      t0 = 0.5*p2
      h0 = ga0 + t0
    echo "Begin H ",nsnow," : ",h0,"  Sg: ",ga0,"  T: ",t0," g.norm2: ",gs[nsnow].norm2

    if n >= qnbegin:
      md.evolve qntau
      hinvp.invH(lbfgs, p)
      threads:
        var p2t = 0.0
        for i in 0..<p.len:
          p2t += p[i].redot hinvp[i]
        threadMaster: p2 = p2t
    else:
      md.evolve tau
      threads:
        var p2t = 0.0
        for i in 0..<p.len:
          p2t += p[i].norm2
        threadMaster: p2 = p2t

    let
      ga1 = gs[nsnow].gaugeAction2 gc
      t1 = 0.5*p2
      h1 = ga1 + t1
    echo "End H ",nsnow," : ",h1,"  Sg: ",ga1,"  T: ",t1," g.norm2: ",gs[nsnow].norm2

    #when true:
    when false:
      block:
        var g1 = lo.newgauge
        var p1 = lo.newgauge
        threads:
          for i in 0..<g1.len:
            g1[i] := gs[nsnow][i]
            p1[i] := p[i]
            p[i] := -1*p[i]
        H.evolve tau
        threads:
          var p2t = 0.0
          for i in 0..<p.len:
            p2t += p[i].norm2
          threadMaster: p2 = p2t
        let
          ga1 = gs[nsnow].gaugeAction2 gc
          t1 = 0.5*p2
          h1 = ga1 + t1
        echo "Reversed H: ",h1,"  Sg: ",ga1,"  T: ",t1
        echo "Reversibility: dH: ",h1-h0,"  dSg: ",ga1-ga0,"  dT: ",t1-t0
        #echo p[0][0]
        for i in 0..<g1.len:
          gs[nsnow][i] := g1[i]
          p[i] := p1[i]

    let
      dH = h1 - h0
      acc = exp(-dH)
      accr = R.uniform
    if accr <= acc:  # accept
      echo "M-H ",nsnow," : ACCEPT  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
    else:  # reject
      echo "M-H ",nsnow," : REJECT  dH: ",dH,"  exp(-dH): ",acc,"  r: ",accr
      threads:
        for i in 0..<gs[nsnow].len:
          gs[nsnow][i] := g0[i]

    let pl = gs[nsnow].plaq3
    echo "plaq ",nsnow," : ",pl.re," ",pl.im
    echo "topo ",nsnow," : ",gs[nsnow].topo2DU1

    if n >= qnbegin:
      if ns < nstream-1:
        f.getforce gs[nsnow]
        lbfgs.add(x = gs[nsnow], f = f)
      else:  # Skip add at the end, only reverse
        f.getforce gs[nstream-1-nsnow]
        lbfgs.reverseadd(x = gs[nstream-1-nsnow], f = f)
        forward = not forward

echoTimers()
checkgc()
qexfinalize()
