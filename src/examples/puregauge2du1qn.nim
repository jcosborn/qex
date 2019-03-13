import qex
import gauge, physics/qcdTypes
import mdevolve
import os, strutils, algorithm, macros, times
import random  # for seeding our RNGs

const CHECK = false
const CHECKLBFGS = false
const CHECKLBFGSMOM = false
const CHECKLBFGSEIGEN = false and not CHECKLBFGSMOM
const CHECKHESSIAN = true
const CHECKREVERSIBLE = false

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

type LBFGS[F] = ref object
  ys,gamma:seq[float]
  y,s,v,u:seq[F]  ## v, u for the product form
  sortedix:seq[int]
  p:int  ## Index for a new addition, or 1 plus the end of the ring.  Temporary storage at p-1.
  sortedlen:int
  lambda:float  ## minimum eigen value of H to regularize the zero modes.
  invsqrtH0:float
  yscale:float  ## scale y <- yscale*y
  first:bool

proc initLBFGS(lo:Layout, n:int, h0 = 1.0, lambda:float = 0.0, yscale:float = 1.0):auto =
  type F = type(lo.newGauge)
  var r {.noinit.} :LBFGS[F]
  r.new
  r.p = 0
  r.first = true
  r.ys.newseq(n)
  r.y.newseq(n)
  r.s.newseq(n)
  r.sortedlen = 0
  r.sortedix.newseq(n-2)  # !IMPORTANT!  N states -> N differences ->  N-2 excluding 1 state
  r.v.newseq(n)
  r.u.newseq(n)
  r.gamma.newseq(n)  # coefficients for inversion, g = v^\dag u - 1
  r.invsqrtH0 = 1.0 / sqrt(h0)
  r.lambda = lambda
  r.yscale = yscale
  for i in 0..<n:
    r.y[i] = lo.newGauge
    r.s[i] = lo.newGauge
  for i in 0..<n:
    r.v[i] = lo.newGauge
    r.u[i] = lo.newGauge
  r

proc add[F](o:LBFGS[F]; x,f:F) =
  var n = o.p
  if o.first:
    o.first = false
  else:
    # save new s = ln(x_old x.adj), y = yscale*(f_old-f)
    let n1 = if n == 0: o.ys.len - 1 else: n - 1
    threads:
      var ys = 0.0
      for mu in 0..<x.len:
        for i in x[mu]:
          when CHECKLBFGS:
            var t = o.s[n1][mu][i] - x[mu][i]
          else:
            var t = ln(o.s[n1][mu][i] * x[mu][i].adj)
            t.projectTAH t
            #t := t.im
          when CHECK:
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
          #when CHECKLBFGS:
          #  var r = o.y[n1][mu][i] - f[mu][i]
          #else:
          #  var r = (o.y[n1][mu][i] - f[mu][i]).im
          #r *= o.yscale
          let r = o.yscale*(o.y[n1][mu][i] - f[mu][i])
          when CHECK:
            if i == 0:
              let
                yi = o.y[n1][mu][i]
                fi = f[mu][i]
              echo "fi_old: ",yi
              echo "fi: ",fi
              echo "fi-fi_old: ",fi-yi
              echo "==r: ",r
          o.y[n1][mu][i] := r
        ys += o.y[n1][mu].redot o.s[n1][mu]
      threadMaster:
        o.ys[n1] = ys
    when CHECK:
      let ysd = o.y[n1].dot o.s[n1]
      echo "ysd: ",ysd
      echo "ys[",n1,"]: ",o.ys[n1].float
  for mu in 0..<x.len:
    o.s[n][mu] := x[mu]  # TODO: we can elide this extra copy
    o.y[n][mu] := f[mu]
    if CHECKLBFGS:
      echo "    x[mu][0]: ",x[mu][0]
      echo "    f[mu][0]: ",f[mu][0]
  inc n
  if n == o.ys.len: n = 0
  o.p = n

proc reverseadd[F](o:LBFGS[F]; x,f:F) =
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
  threads:
    for i in 0..<(p div 2):
      for mu in 0..<o.y[p].len:  # FIXME: use one more index indirection to avoid copy
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

proc A[F](o:LBFGS[F], k:int, z:F) =
  ## A_k z  ->  z  where  H_k = A_k A_k^\dag, A_k = (1 - u v^\dag) A_{k-1}
  threads:
    # A0 z -> z
    let a0 = 1.0 / o.invsqrtH0
    for mu in 0..<z.len: z[mu] *= a0
    for j in 0..k:
      let i = o.sortedix[j]
      let vz = -(o.v[i].dot z)
      for mu in 0..<z.len:
        for e in z[mu]:
          let t = z[mu][e] + vz * o.u[i][mu][e]
          z[mu][e] := t
proc Adag[F](o:LBFGS[F], k:int, z:F) =
  ## A_k^\dag z  ->  z  where  H_k = A_k A_k^\dag, A_k^\dag = A_{k-1}^\dag (1 - v u^\dag)
  threads:
    for j in countdown(k,0):
      let i = o.sortedix[j]
      let uz = -(o.u[i].dot z)
      for mu in 0..<z.len:
        for e in z[mu]:
          let t = z[mu][e] + uz * o.v[i][mu][e]
          z[mu][e] := t
    # A0d z -> z
    let a0 = 1.0 / o.invsqrtH0
    for mu in 0..<z.len: z[mu] *= a0
proc H[F](o:LBFGS[F], k:int, r:F) =
  when CHECKLBFGS and false: echo "r: ",r
  o.Adag(k, r)
  when CHECKLBFGS and false: echo "Adr: ",r
  o.A(k, r)
  when CHECKLBFGS and false: echo "AAdr: ",r
proc H[F](r:F, o:LBFGS[F], k:int, z:F) =
  for mu in 0..<z.len: r[mu] := z[mu]
  o.H(k,r)
proc H[F](o:LBFGS[F], r:F) =
  let n = o.sortedlen
  o.H(n-1, r)
proc H[F](r:F, o:LBFGS[F], z:F) =
  let n = o.sortedlen
  r.H(o, n-1, z)
proc sqrtH[F](o:LBFGS[F]; x:F) =
  ## A x  ->  x  where  H = A A^\dag
  let n = o.sortedlen
  o.A(n-1, x)

proc B[F](o:LBFGS[F], k:int, z:F) =
  ## B_k z  ->  z  where  H_k^-1 = B_k B_k^\dag, B_k = (1 - v u^\dag / g) B_{k-1}
  threads:
    # B0 z -> z
    let b0 = o.invsqrtH0
    for mu in 0..<z.len: z[mu] *= b0
    for j in 0..k:
      let i = o.sortedix[j]
      let uz = (o.u[i].dot z) / (-o.gamma[i])
      for mu in 0..<z.len:
        for e in z[mu]:
          let t = z[mu][e] + uz * o.v[i][mu][e]
          z[mu][e] := t
proc Bdag[F](o:LBFGS[F], k:int, z:F) =
  ## A_k^\dag z  ->  z  where  H_k^-1 = B_k B_k^\dag, B_k^\dag = B_{k-1}^\dag (1 - u v^\dag / g)
  threads:
    for j in countdown(k,0):
      let i = o.sortedix[j]
      let vz = (o.v[i].dot z) / (-o.gamma[i])
      for mu in 0..<z.len:
        for e in z[mu]:
          let t = z[mu][e] + vz * o.u[i][mu][e]
          z[mu][e] := t
    # B0d z -> z
    let b0 = o.invsqrtH0
    for mu in 0..<z.len: z[mu] *= b0
proc invH[F](o:LBFGS[F], k:int, r:F) =
  when CHECKLBFGS and false: echo "r: ",r
  o.Bdag(k, r)
  when CHECKLBFGS and false: echo "Bdr: ",r
  o.B(k, r)
  when CHECKLBFGS and false: echo "BBdr: ",r
proc invH[F](r:F, o:LBFGS[F], k:int, z:F) =
  for mu in 0..<z.len: r[mu] := z[mu]
  o.invH(k,r)
proc invH[F](o:LBFGS[F], r:F) =
  let n = o.sortedlen
  o.invH(n-1, r)
proc invH[F](r:F, o:LBFGS[F], z:F) =
  let n = o.sortedlen
  r.invH(o, n-1, z)
proc sqrtInvH[F](o:LBFGS[F]; x:F) =
  ## B x  ->  x  where  H^-1 = B B^\dag
  let n = o.sortedlen
  o.B(n-1, x)

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
        threads: g.fromPrimmeArrayGauge xp
        m.lbfgs.op g
        threads: g.toPrimmeArrayGauge yp
      err[] = 0
  genMatvec(invH)
  genMatvec(H)
  genMatvec(sqrtH)

proc prep[F](o:LBFGS[F], cutoff = 0.0, reduce = 0) =
  ## For differences y and s, the item k has the difference between k and k+1.
  ## Before using the approximation, o.p points to the current updating stream number.
  ## We sort according to ys, excluding item o.p and o.p-1, so we don't depend on ourselves.
  let t0 = epochTime()
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
  for j in 0..<o.sortedix.len:
    let i = o.sortedix[j]
    var yy,ss:float
    threads:
      let y2 = o.y[i].norm2
      let s2 = o.s[i].norm2
      threadMaster:
        yy = y2
        ss = s2
    when CHECK:
      echo j," ix ",i," ys ",o.ys[i]," ys/yy ",(o.ys[i]/yy)," yy ",yy," ss ",ss
    if unset and o.ys[i] <= cutoff:
      o.sortedlen = j
      unset = false
  if reduce > 0:
    echo "lBGFS nmem reduce: ",reduce
    o.sortedlen -= reduce
    if o.sortedlen < 0: o.sortedlen = 0
  echo "lBFGS nmem = ",o.sortedlen
  for j in 0..<o.sortedlen:  # Compute v_i, u_i, and gamma_i
    let i = o.sortedix[j]
    o.u[i].H(o, j-1, o.s[i])  # temporarily u <- G_{k-1} s_k  (j is the actual order number)
    o.v[i].invH(o, j-1, o.y[i])  # temporarily v <- G_{k-1}^{-1} y_k
    var ss,sgs,ygiy:float
    threads:
      let
        tss = o.s[i].norm2
        tsgs = o.s[i].redot o.u[i]
        tygiy = o.y[i].redot o.v[i]
      threadMaster:
        ss = tss
        sgs = tsgs
        ygiy = tygiy
    var delta1 = o.lambda*ss/sgs  # = 1-delta
    if delta1 > 1.0: delta1 = 1.0
    let
      delta = 1.0 - delta1
      wgiw = ygiy / o.ys[i]
      cy = 1.0 / sqrt(o.ys[i])
      cs = sqrt(delta/sgs)
      wgiz = cs/cy
    o.gamma[i] = sqrt(delta1*(wgiw - o.ys[i]/sgs + 1.0) + o.ys[i]/sgs)  # negative sign works, too
    let
      theta = (delta1+o.gamma[i]-wgiz) / (2.0*wgiz+wgiw+delta)
      cyv = cy*theta
      csv = cs*(1.0+theta)
    when CHECK:
      echo j," delta: ",delta," sgs: ",sgs," ygiy: ",ygiy," gamma: ",o.gamma[i].float
    threads:
      for mu in 0..<o.u[i].len:
        for e in o.u[i][mu]:
          let t = cy*o.y[i][mu][e] + cs*o.u[i][mu][e]
          o.u[i][mu][e] := t
      for mu in 0..<o.v[i].len:
        for e in o.v[i][mu]:
          let t = cyv*o.v[i][mu][e] + csv*o.s[i][mu][e]
          o.v[i][mu][e] := t
  let t1 = epochTime()
  echo "LBFGS prep time: ",t1-t0
  when CHECKLBFGSEIGEN:
    let te0 = epochTime()
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
        vecs = newAlignedMemU[ccomplex[float]]int(p.numEvals*p.nLocal)
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
        let t = vecs[0]  # XXX Force instantiating the type, works around Nim compiler bug.
        let ret = p.run(vals, vecs[0].addr, rnorms)
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
        echo i," ev ",vals[i].float," rnorm ",rnorms[i].float
    o.lbfgseigen(primme_smallest)
    o.lbfgseigen(primme_largest)
    let te1 = epochTime()
    echo "LBFGS eigen time: ",te1-te0
    #qexExit 0

when CHECKLBFGS:
  const cklbfgsnd = 32
  let nseq = 8192
  type V = array[cklbfgsnd,float]
  type rV = ref[V]
  type M = array[cklbfgsnd,V]
  type D = distinct rV
  proc initD:D = D(V.new)
  template newOneOf(x:D):D = initD()
  template len(x:D):int = x.rV[].len
  template `[]`(x:D,i:untyped):untyped = x.rV[][i]
  template `[]=`(x:D,i:untyped,y:untyped):untyped = x.rV[][i] = y
  template `$`(x:D):string = $x.rV[]
  iterator items(x:D):int =
    for i in 0..<x.len: yield i
  template getdiff(z,x,y:D) =
    for i in x:
      z[i] = x[i] - y[i]
  proc dot(x,y:D):float =
    for i in x:
      result += x[i]*y[i]
  template `*=`(x:D,y:float) =
    for i in x:
      x[i] *= y
  template `-=`(x:D,y:D) =
    for i in x:
      x[i] -= y[i]
  template `:=`(x,y:float) = x = y
  template `:=`(x,y:D) =
    for i in x:
      x[i] = y[i]
  template redot(x,y:D):float = x.dot y
  template norm2(x:D):float = x.dot x
  type F = seq[D]
  proc dot(x,y:F):float =
    for i in 0..<x.len:
      result += x[i].dot y[i]
  template redot(x,y:F):float = x.dot y
  template norm2(y:F):float =
    let x = y
    x.dot x
  template getforce(z:D,m:M,y:D)=
    for j in 0..<m.len:
      for i in 0..<m.len:
        z[i] += m[i][j] * y[j]
  proc checklbfgs =
    var cklbfgsr {.noinit.} :LBFGS[F]
    cklbfgsr.new
    cklbfgsr.p = 0
    cklbfgsr.first = true
    cklbfgsr.ys.newseq(nseq)
    cklbfgsr.y.newseq(nseq)
    cklbfgsr.s.newseq(nseq)
    cklbfgsr.sortedlen = 0
    cklbfgsr.sortedix.newseq(nseq-2)  # N states -> N differences ->  N-2 excluding 1 state
    cklbfgsr.v.newseq(nseq)
    cklbfgsr.u.newseq(nseq)
    cklbfgsr.gamma.newseq(nseq)
    cklbfgsr.invsqrtH0 = 1.0
    cklbfgsr.yscale = 1.0
    cklbfgsr.lambda = 0.0
    for i in 0..<nseq:
      cklbfgsr.y[i].newseq(1)
      cklbfgsr.y[i][0] = initD()
      cklbfgsr.s[i].newseq(1)
      cklbfgsr.s[i][0] = initD()
    for i in 0..<nseq:
      cklbfgsr.v[i].newseq(1)
      cklbfgsr.v[i][0] = initD()
      cklbfgsr.u[i].newseq(1)
      cklbfgsr.u[i][0] = initD()
    var rng = initRand 5^5
    var m:M
    for i in 0..<m.len-1:
      m[i][i] = float((i+1)*(i+1))
      m[i][i+1] = 1
      m[i+1][i] = 1
    m[m.len-1][m.len-1] = 1
    for i in 0..<nseq:
      var x = initD()
      for k in 0..<cklbfgsnd: x[k] = rng.rand 1.0
      #echo "x: ",x
      var y = initD()
      y.getforce(m, x)
      #echo "y: ",y
      cklbfgsr.add(x = @[x], f = @[y])
    cklbfgsr.prep
    var x = initD()
    for i in x: x[i] = 1
    echo "test x: ",x
    var y = initD()
    y.getforce(m, x)
    echo "test y: ",y
    var hx = @[x]
    var tmp = initD()
    cklbfgsr.H hx
    tmp.getdiff(hx[0], y)
    echo "hx: ",hx
    echo "hx err: ",norm2(tmp)
    var hinvhx = hx
    cklbfgsr.invH hinvhx
    tmp.getdiff(hinvhx[0], x)
    echo "hinvhx: ",hinvhx
    echo "hinvhx err: ",norm2(tmp)
    var hinvy = @[y]
    cklbfgsr.invH hinvy
    tmp.getdiff(hinvy[0], x)
    echo "hinvy: ",hinvy
    echo "hinvy err: ",norm2(tmp)
    var hsqrtx = @[x]
    cklbfgsr.sqrtH hsqrtx
    echo "hsqrtx: ",hsqrtx
  qexinit()
  checklbfgs()
  qexExit 0

qexinit()
threads: echo "thread ",threadNum," / ",numThreads

let
  beta = floatParam("beta", 5.0)
  nx = intParam("nx", 64)
  nt = intParam("nt", nx)
  trajs = intParam("trajs", 512)
  qnbegin = intParam("qnbegin",64)
  tau = floatParam("tau", 2.0)
  steps = intParam("steps", 10)
  qntau = floatParam("qntau", tau)
  qnsteps = intParam("qnsteps", steps)
  qnyscut = floatParam("qnyscut",0)
  qnyscale = floatParam("qnyscale", 1.0/(2.0*beta))  # scaling the approximated Hessian by inverse diagonal term of the free field Hessian
  qnh0 = floatParam("qnh0",1.0)  # initial diagonal term of the approximated Hessian
  seed = intParam("seed", 11^11)
  gfix = intParam("gfix", 1).bool
  gfixextra = intParam("gfixextra", 0).bool
  gfixunit = intParam("gfixunit", 1).bool
  nstream = intParam("nstream", 10)
  lambda = floatParam("lambda", 0.1)
  randomInit = intParam("randomInit", 1).bool

echo "beta = ",beta
echo "nx = ",nx
echo "nt = ",nt
echo "trajs = ",trajs
echo "tau = ",tau
echo "steps = ",steps
echo "qnbegin = ",qnbegin
echo "qntau = ",qntau
echo "qnsteps = ",qnsteps
echo "qnyscut = ",qnyscut
echo "qnyscale = ",qnyscale
echo "qnh0 = ",qnh0
echo "gfix = ",gfix.int
echo "gfixextra = ",gfixextra.int
echo "gfixunit = ",gfixunit.int
echo "nstream = ",nstream
echo "lambda = ",lambda
echo "seed = ",seed
echo "randomInit = ",randomInit.int

var (gs,r,R) = setup([nx,nt],nstream,seed)

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
  lbfgs = lo.initLBFGS(nstream, h0 = qnh0, lambda = lambda, yscale = qnyscale)
  p = lo.newgauge
  f = lo.newgauge
  g0 = lo.newgauge
  nsnow = 0  # The serial number of the current running HMC.
  forward = true
  hinvp = lo.newgauge

for i in 0..<nstream:
  if randomInit:
    gs[i].random r
  else:
    for mu in 0..<gs[i].len: gs[i][mu] := 1
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

when CHECKLBFGSMOM:
  proc checkmom =
    # Refreshing the momentum and check lbfgs stability.
    proc checklbfgsmom =
      var p2 = 0.0
      threads: p.randomTAH r
      if gfix: p.maxTreeFix(0.0, gfixextra)
      threads:
        var p2t = 0.0
        for i in 0..<p.len:
          p2t += p[i].norm2
        threadMaster: p2 = p2t
      lbfgs.sqrtH p
      #if gfix: p.maxTreeFix(0.0, gfixextra)  # FIXME
      hinvp.invH(lbfgs, p)
      #threads:
      block:
        var p2t,p2nt,hp2 = 0.0
        for i in 0..<p.len:
          p2t += p[i].redot hinvp[i]
          p2nt += hinvp[i].norm2
          hp2 += p[i].norm2
        #threadMaster:
        block:
          let e = abs(p2-p2t)/max(p2.abs,p2t.abs)
          echo "SqrtHP.norm2: ",hp2," HinvP.norm2: ",p2nt
          echo "p2: ",p2
          echo "p2_Hinv: ",p2t," err: ",e
    for c in countdown(10,0):
      let cut = float(c*lat[0])
      echo "LBFGS cutoff: ",cut
      lbfgs.prep(cutoff = cut)
      for i in 0..<1024: checklbfgsmom()

when CHECKHESSIAN:
  proc checkhession(nsnow:int) =
    # Goodness of the Hessian approximation
    # ./bin/puregauge2du1qn -beta:32 -nx:8 -steps:20 -gfix:1 -gfixextra:1
    #     -nstream:400 -qnbegin:8 -trajs:9 -randomInit:0 -lambda:0
    # Should see ratios of errors with and without the 2nd order term about -0.0204(8)
    let dev = 0.01
    let s0 = gs[nsnow].gaugeAction2 gc
    f.getforce gs[nsnow]
    for i in 0..<g0.len:
      for j in g0[i]:
        let t = dev * ln gs[nsnow][i][j]
        g0[i][j] := t
    if gfix: g0.maxTreeFix(0.0, gfixextra)
    let s01 = f.redot g0
    p.H(lbfgs, g0)
    let s02 = (0.5/lbfgs.yscale) * p.redot g0
    let s1a = s0 + s01 + s02
    for i in 0..<g0.len:
      for j in g0[i]:
        let t = gs[nsnow][i][j] * exp(dev * ln gs[nsnow][i][j])
        g0[i][j] := t
    if gfix: g0.maxTreeFix(1.0, gfixextra)
    let s1 = g0.gaugeAction2 gc
    echo "s0: ",s0," s01: ",s01," s02: ",s02
    echo "s1a: ",s1a
    echo "s1: ",s1
    let e12 = s1-s1a
    let e11 = s1-(s0+s01)
    echo "error: ",e12, " s1-(s0+s01): ",e11, " ratio: ",e12/e11

proc getp2:float =
  var p2 = 0.0
  threads:
    var p2t = 0.0
    for i in 0..<p.len:
      p2t += p[i].norm2
    threadMaster:
      p2 = p2t
  p2

proc getpgp:float =
  var p2 = 0.0
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
  p2

when CHECKREVERSIBLE:
  proc checkreversible(nsnow:int) =
    var g1 = lo.newgauge
    var p1 = lo.newgauge
    threads:
      for i in 0..<g1.len:
        g1[i] := gs[nsnow][i]
        p1[i] := p[i]
        p[i] := -1*p[i]
    H.evolve tau
    p2 = getp2()
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

for n in 1..trajs:
  let tt0 = epochTime()
  if n == qnbegin:
    echo "STARTING QN update"
    md.steps = qnsteps
    for i in 0..<nstream:
      f.getforce gs[i]
      lbfgs.add(x = gs[i], f = f)
    when CHECKLBFGSMOM:
      checkmom()
      echoTimers()
      checkgc()
      qexExit 0

  for ns in 0..<nstream:
    let ts0 = epochTime()
    nsnow = if forward: ns else: nstream-1-ns
    echo "Begin traj: ",n," nsNow: ",nsnow," ",(if forward: "forward" else: "backward")
    if n >= qnbegin:
      lbfgs.prep(cutoff = qnyscut * lat[0].float)
      when CHECKHESSIAN:
        checkhession(nsnow)
        #qexexit()

    var p2 = 0.0
    threads: p.randomTAH r
    if gfix: p.maxTreeFix(0.0, gfixextra)
    if CHECK or n < qnbegin:
      p2 = getp2()

    if n >= qnbegin:
      lbfgs.sqrtH p
      hinvp.invH(lbfgs, p)
      p2 = getpgp()

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
      p2 = getpgp()
    else:
      md.evolve tau
      p2 = getp2()

    let
      ga1 = gs[nsnow].gaugeAction2 gc
      t1 = 0.5*p2
      h1 = ga1 + t1
    echo "End H ",nsnow," : ",h1,"  Sg: ",ga1,"  T: ",t1," g.norm2: ",gs[nsnow].norm2

    when CHECKREVERSIBLE:
      checkreversible(nsnow)

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

    let ts1 = epochTime()
    echo "Stream trajectory time: ",ts1-ts0

  let tt1 = epochTime()
  echo "Total trajectory time: ",tt1-tt0

echoTimers()
checkgc()
qexfinalize()
