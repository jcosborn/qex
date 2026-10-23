import base
import layout
import layout/shifts
#import profile
#import stdUtils
#import gaugeUtils

proc newLinkShifts*[V:static[int],T](lo: Layout[V], t: typedesc[T], len: int): seq[seq[ShiftB[T]]] =
  ## s[mu][nu], mu != nu: a shift along nu by len for fields of lo with site type T.
  ## A ShiftB holds indices and buffers only, so a row serves any such field.
  var l = lo  # getShift caches in the layout
  let nd = l.nDim
  result.newSeq(nd)
  for mu in 0..<nd:
    result[mu].newSeq(nd)
    for nu in 0..<nd:
      if mu != nu: result[mu][nu].initShiftB(l, T, nu, len, "all")

proc startCornerShifts*[V:static[int],T](uu: openArray[Field[V,T]]): auto =
  ## s[mu][nu] <- U_mu(x+nu) for mu != nu. This call and the kernel completing
  ## s read uu through a raw pointer while allocating; the caller keeps uu and
  ## its fields referenced until that kernel returns.
  let u = cast[ptr cArray[Field[V,T]]](unsafeAddr(uu[0]))
  let nd = uu.len
  var s = newLinkShifts(uu[0].l, T, 1)
  threads:
    for mu in 0..<nd:
      for nu in 0..<nd:
        if mu!=nu:
          s[mu][nu].startSB(u[mu][ix])
  s

proc startPlaqShifts*[T](u0: openArray[T]): auto =
  var u = @u0
  var sf:seq[seq[ShiftB[type(u[0][0])]]]
  var sb:seq[seq[ShiftB[type(u[0][0])]]]
  let nd = u.len
  sf.newSeq(nd)
  sb.newSeq(nd)
  for mu in 0..<nd:
    sf[mu].newSeq(nd)
    sb[mu].newSeq(nd)
    for nu in 0..<nd:
      sb[mu][nu].initShiftB(u[mu], nu, -1, "all")
      if mu!=nu:
        sf[mu][nu].initShiftB(u[mu], nu, 1, "all")
  threads:
    for mu in 0..<nd:
      sb[mu][mu].startSB(u[mu][ix])
    for mu in 0..<nd:
      for nu in 0..<nd:
        if mu!=nu:
          sf[mu][nu].startSB(u[mu][ix])
          sb[mu][nu].startSB(u[mu][ix])
  return (sf,sb)

proc finishPlaqShifts*[T](u0,b0: openArray[T], sf,sb: seq): auto =
  var u = @u0
  var b = @b0
  var sc:seq[seq[ShiftB[type(u[0][0])]]]
  let nd = u.len
  sc.newSeq(nd)
  for mu in 0..<nd:
    sc[mu].newSeq(nd)
    for nu in 0..<nd:
      if mu!=nu:
        sc[mu][nu].initShiftB(u[mu], nu, 1, "all")
  threads:
    for mu in 0..<nd:
      for ir in b[mu]:
        localSB(sb[mu][mu], ir, assign(b[mu][ir], it), u[mu][ix])
    for mu in 0..<nd:
      boundarySB(sb[mu][mu], assign(b[mu][ir], it))
    threadBarrier()
    for mu in 0..<nd:
      for nu in 0..<nd:
        if mu!=nu:
          sc[mu][nu].startSB(b[mu][ix])
  return (sc)

proc startStapleShifts*[T](u: openArray[T]): auto =
  var s:seq[seq[seq[ShiftB[type(u[0][0][0])]]]]
  let nd = u.len
  s.newSeq(nd)
  for mu in 0..<nd:
    s[mu].newSeq(nd)
    for nu in 0..<nd:
      if nu!=mu:
        s[mu][nu].newSeq(nd)
        for sig in 0..<nd:
          if sig!=mu:
            s[mu][nu][sig].initShiftB(u[mu][nu], sig, 1, "all")
            s[mu][nu][sig].startSB(u[mu][nu][ix])
  return s

proc stapleFields[V:static[int],T](u: openArray[Field[V,T]], reuse: FieldArray[V,T]): auto =
  let lo = u[0].l
  let nd = lo.nDim
  if reuse.shape.len == 0:
    return newFieldArray2(lo, type(u[0]), [nd,nd], mu != nu)
  if reuse.shape.len != 2 or reuse.shape[0] != nd or reuse.shape[1] != nd:
    raise newException(ValueError, "staple storage requires its original shape")
  for i in 0..<nd*nd:
    let f = reuse[i]
    if f != nil:
      if f.l != lo:
        raise newException(ValueError, "staple storage requires its original layout")
      for x in u:
        if f.s.mem != nil and f.s.mem == x.s.mem:
          raise newException(ValueError, "staple storage must be disjoint from inputs")
  reuse

# default(FieldArray[V,T]) arguments break Nim 2.0; use the method form.
proc makeFwdStaples*[V:static[int],T](uu: openArray[Field[V,T]], s: auto,
    reuse: FieldArray[V,T] = FieldArray[V,T].default): auto =
  ## st[mu,nu](x) = U_nu(x) U_mu(x+nu) U_nu(x+mu)†.
  ## reuse: an array this proc returned for uu's layout, with every mu != nu entry
  ## allocated and storage disjoint from uu. Shape, layout and aliasing raise
  ## ValueError; a missing entry crashes in the kernel.
  mixin mul, load1, adj
  tic()
  let u = cast[ptr cArray[type(uu[0])]](unsafeAddr(uu[0]))
  let lo = u[0].l
  let nd = lo.nDim
  let nc = u[0][0].ncols
  let flops = lo.nSites.float*float(nd*(nd-1)*3*(4*nc-1)*nc*nc)
  var st = stapleFields(uu, reuse)
  toc("makeStaples setup")
  threads:
    tic()
    var umu,unu,umunu: type(load1(u[0][0]))
    for ir in lo:
      for mu in 1..<nd:
        for nu in 0..<mu:
          if isLocal(s[mu][nu],ir) and isLocal(s[nu][mu],ir):
            localSB(s[mu][nu], ir, assign(umu,it), u[mu][ix])
            localSB(s[nu][mu], ir, assign(unu,it), u[nu][ix])
            mul(umunu, umu, unu.adj)
            mul(st[mu,nu][ir], u[nu][ir], umunu)
            mul(st[nu,mu][ir], u[mu][ir], umunu.adj)
    toc("makeStaples local")
    #[
    var needBoundary = false
    for mu in 0..<nd:
      for nu in 0..<nd:
        if mu != nu:
          boundaryWaitSB(s[mu][nu]): needBoundary = true
    toc("makeStaples wait")
    if needBoundary:
      boundarySyncSB()
      for ir in lo:
        for mu in 1..<nd:
          for nu in 0..<mu:
            if not isLocal(s[mu][nu],ir) or not isLocal(s[nu][mu],ir):
              getSB(s[mu][nu], ir, assign(umu,it), u[mu][ix])
              getSB(s[nu][mu], ir, assign(unu,it), u[nu][ix])
              mul(umunu, umu, unu.adj)
              mul(st[mu,nu][ir], u[nu][ir], umunu)
              mul(st[nu,mu][ir], u[mu][ir], umunu.adj)
    ]#
    for mu in 1..<nd:
      for nu in 0..<mu:
        var needBoundary = false
        boundaryWaitSB(s[mu][nu]): needBoundary = true
        boundaryWaitSB(s[nu][mu]): needBoundary = true
        if needBoundary:
          # Thread 0 waits for rbuf; synchronize before getSB.
          boundarySyncSB()
          for ir in lo:
            if not isLocal(s[mu][nu],ir) or not isLocal(s[nu][mu],ir):
              getSB(s[mu][nu], ir, assign(umu,it), u[mu][ix])
              getSB(s[nu][mu], ir, assign(unu,it), u[nu][ix])
              mul(umunu, umu, unu.adj)
              mul(st[mu,nu][ir], u[nu][ir], umunu)
              mul(st[nu,mu][ir], u[mu][ir], umunu.adj)
    toc("makeStaples boundary")
  toc("makeStaples threads", flops=flops)
  return st

proc staplePlanes*(nd: int): seq[(int,int)] =
  ## The planes (mu, nu), nu < mu, in lexicographic order.
  for mu in 1..<nd:
    for nu in 0..<mu: result.add (mu,nu)

proc makeStaples*[V:static[int],T](uu: openArray[Field[V,T]], s: auto,
    fwd: FieldArray[V,T] = FieldArray[V,T].default,
    bck: FieldArray[V,T] = FieldArray[V,T].default; post = true): auto =
  ## sft: fwd staples
  ## stu: bck staples, offset up
  ## ss: shifts for stu, started when post; the caller completes them.
  ## fwd, bck: arrays this proc returned for uu's layout, with every mu != nu entry
  ## allocated and storage disjoint from uu and from each other. Shape, layout and
  ## aliasing raise ValueError; a missing entry crashes in the kernel.
  mixin mul, load1, adj
  tic()
  let u = cast[ptr cArray[type(uu[0])]](unsafeAddr(uu[0]))
  let lo = u[0].l
  let nd = lo.nDim
  let nc = u[0][0].ncols
  let flops = lo.nSites.float*float(nd*(nd-1)*6*(4*nc-1)*nc*nc)
  var
    stf = stapleFields(uu, fwd)
    stu = stapleFields(uu, bck)
  if nd > 1 and stf[0,1].s.mem == stu[0,1].s.mem:
    raise newException(ValueError, "forward and backward staples must be disjoint")
  var ss = newLinkShifts(lo, T, -1)
  toc("makeStaples setup")
  threads:
    tic()
    var umu,unu,umunu,unumu {.noInit.}: evalType(u[0][0])
    for ir in lo:
      for mu in 1..<nd:
        for nu in 0..<mu:
          if isLocal(s[mu][nu],ir) and isLocal(s[nu][mu],ir):
            localSB(s[mu][nu], ir, assign(umu,it), u[mu][ix])
            localSB(s[nu][mu], ir, assign(unu,it), u[nu][ix])
            mul(umunu, umu, unu.adj)
            mul(stf[mu,nu][ir], u[nu][ir], umunu)
            mul(stf[nu,mu][ir], u[mu][ir], umunu.adj)
            mul(unumu, u[nu][ir].adj, u[mu][ir])
            mul(stu[mu,nu][ir], unumu, unu)
            mul(stu[nu,mu][ir], unumu.adj, umu)
    toc("makeStaples local")
    var needBoundary = false
    for mu in 1..<nd:
      for nu in 0..<mu:
        var needBoundaryU = false
        boundaryWaitSB(s[mu][nu]): needBoundaryU = true
        boundaryWaitSB(s[nu][mu]): needBoundaryU = true
        needBoundary = needBoundary or needBoundaryU
        if needBoundaryU:
          boundarySyncSB()
          for ir in lo:
            if not isLocal(s[mu][nu],ir) or not isLocal(s[nu][mu],ir):
              getSB(s[mu][nu], ir, assign(umu,it), u[mu][ix])
              getSB(s[nu][mu], ir, assign(unu,it), u[nu][ix])
              mul(unumu, u[nu][ir].adj, u[mu][ir])
              mul(stu[mu,nu][ir], unumu, unu)
              mul(stu[nu,mu][ir], unumu.adj, umu)
        if post:
          threadBarrier()
          ss[mu][nu].startSB(stu[mu,nu][ix])
          ss[nu][mu].startSB(stu[nu,mu][ix])
    toc("makeStaplesU boundary")
    if needBoundary:
      boundarySyncSB()
      for ir in lo:
        for mu in 1..<nd:
          for nu in 0..<mu:
            if not isLocal(s[mu][nu],ir) or not isLocal(s[nu][mu],ir):
              getSB(s[mu][nu], ir, assign(umu,it), u[mu][ix])
              getSB(s[nu][mu], ir, assign(unu,it), u[nu][ix])
              mul(umunu, umu, unu.adj)
              mul(stf[mu,nu][ir], u[nu][ir], umunu)
              mul(stf[nu,mu][ir], u[mu][ir], umunu.adj)
    toc("makeStaplesF boundary")
  toc("makeStaples threads", flops=flops)
  return (stf,stu,ss)


type StapleBatchFields[V:static[int],T] = object
  ## Direction-indexed views of the current scratch batch. Other entries
  ## may alias reused storage and must not be read outside their batch.
  nd: int
  fields: seq[Field[V,T]]

template `[]`(p: StapleBatchFields; mu, nu: int): untyped =
  p.fields[mu*p.nd+nu]

# Half of the six planes in four dimensions: reuse 12 fields for an action or
# a plaquette force and 18 for a force with a rect term, balancing field
# traffic against the cost of fresh allocations.
const temporaryStaplePairs = 3

proc addStaplePlanes*[V:static[int],T](uu: openArray[Field[V,T]], output: auto;
    cp, cr: float; stf, stu, ru: auto;
    cs, ss, sf, sb: seq[seq[ShiftB[T]]]; planes: auto; overwrite = false) =
  ## Shared derivative kernel for full retained staples and one temporary pair.
  ## Planes follow global lexicographic order. overwrite initializes each
  ## direction on its first plane; rectangle scratch is zeroed each batch.
  mixin load1, adj
  let u = cast[ptr cArray[Field[V,T]]](unsafeAddr(uu[0]))
  let f = cast[ptr cArray[type(output[0])]](unsafeAddr(output[0]))
  let lo = u[0].l
  threads:
    tic("staple derivative")
    if cr!=0:
      for plane in planes:
        let (mu, nu) = plane
        sf[mu][nu].startSB(stf[mu,nu][ix])
        sf[nu][mu].startSB(stf[nu,mu][ix])
    if cr!=0:
      for plane in planes:
        let (mu, nu) = plane
        ru[mu,nu] := 0
        ru[nu,mu] := 0
    for ir in u[0]:
      for plane in planes:
        let (mu, nu) = plane
        # plaq
        # In lexicographic plane order, (mu,0) first reaches direction mu;
        # (1,0) first reaches direction zero. Initialize on that contribution.
        if overwrite and nu == 0:
          f[mu][ir] := cp * stf[mu,nu][ir]
        else:
          f[mu][ir] += cp * stf[mu,nu][ir]
        if overwrite and mu == 1:
          f[nu][ir] := cp * stf[nu,mu][ir]
        else:
          f[nu][ir] += cp * stf[nu,mu][ir]
        if isLocal(ss[mu][nu],ir):
          var bmu: type(load1(u[0][0]))
          localSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
          f[mu][ir] += cp * bmu
          if cr!=0:
            var umu,unu,bmunu: type(load1(u[0][0]))
            getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
            getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
            bmunu := bmu * unu
            f[nu][ir] += cr * bmunu * umu.adj
            ru[nu,mu][ir] += bmu.adj * u[nu][ir] * umu
            ru[mu,nu][ir] += u[nu][ir].adj * bmunu
        if isLocal(ss[nu][mu],ir):
          var bnu: type(load1(u[0][0]))
          localSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
          f[nu][ir] += cp * bnu
          if cr!=0:
            var unu,umu,bnumu: type(load1(u[0][0]))
            getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
            getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
            bnumu := bnu * umu
            f[mu][ir] += cr * bnumu * unu.adj
            ru[mu,nu][ir] += bnu.adj * u[mu][ir] * unu
            ru[nu,mu][ir] += u[mu][ir].adj * bnumu
        if cr!=0:
          if isLocal(sf[mu][nu],ir):
            var smu,unu,smunu: type(load1(u[0][0]))
            localSB(sf[mu][nu], ir, assign(smu,it), stf[mu,nu][ix])
            getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
            smunu := smu * unu.adj
            f[mu][ir] += cr * u[nu][ir] * smunu
            f[nu][ir] += cr * u[mu][ir] * smunu.adj
            ru[nu,mu][ir] += u[mu][ir].adj * u[nu][ir] * smu
          if isLocal(sf[nu][mu],ir):
            var snu,umu,snumu: type(load1(u[0][0]))
            localSB(sf[nu][mu], ir, assign(snu,it), stf[nu,mu][ix])
            getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
            snumu := snu * umu.adj
            f[nu][ir] += cr * u[mu][ir] * snumu
            f[mu][ir] += cr * u[nu][ir] * snumu.adj
            ru[mu,nu][ir] += u[nu][ir].adj * u[mu][ir] * snu
    toc("local")
    for plane in planes:
      let (mu, nu) = plane
      var needBoundary = false
      boundaryWaitSB(ss[mu][nu]): needBoundary = true
      boundaryWaitSB(ss[nu][mu]): needBoundary = true
      if needBoundary:
        boundarySyncSB()
        for ir in lo:
          if not isLocal(ss[mu][nu],ir):
            var bmu: type(load1(u[0][0]))
            getSB(ss[mu][nu], ir, assign(bmu,it), stu[mu,nu][ix])
            f[mu][ir] += cp * bmu
            if cr!=0:
              var umu,unu,bmunu: type(load1(u[0][0]))
              getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
              getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
              bmunu := bmu * unu
              f[nu][ir] += cr * bmunu * umu.adj
              ru[nu,mu][ir] += bmu.adj * u[nu][ir] * umu
              ru[mu,nu][ir] += u[nu][ir].adj * bmunu
          if not isLocal(ss[nu][mu],ir):
            var bnu: type(load1(u[0][0]))
            getSB(ss[nu][mu], ir, assign(bnu,it), stu[nu,mu][ix])
            f[nu][ir] += cp * bnu
            if cr!=0:
              var unu,umu,bnumu: type(load1(u[0][0]))
              getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
              getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
              bnumu := bnu * umu
              f[mu][ir] += cr * bnumu * unu.adj
              ru[mu,nu][ir] += bnu.adj * u[mu][ir] * unu
              ru[nu,mu][ir] += u[mu][ir].adj * bnumu
    if cr!=0:
      for plane in planes:
        let (mu, nu) = plane
        var needBoundary = false
        boundaryWaitSB(sf[mu][nu]): needBoundary = true
        boundaryWaitSB(sf[nu][mu]): needBoundary = true
        if needBoundary:
          boundarySyncSB()
          for ir in lo:
            if not isLocal(sf[mu][nu],ir):
              var smu,unu,smunu: type(load1(u[0][0]))
              getSB(sf[mu][nu], ir, assign(smu,it), stf[mu,nu][ix])
              getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
              smunu := smu * unu.adj
              f[mu][ir] += cr * u[nu][ir] * smunu
              f[nu][ir] += cr * u[mu][ir] * smunu.adj
              ru[nu,mu][ir] += u[mu][ir].adj * u[nu][ir] * smu
            if not isLocal(sf[nu][mu],ir):
              var snu,umu,snumu: type(load1(u[0][0]))
              getSB(sf[nu][mu], ir, assign(snu,it), stf[nu,mu][ix])
              getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
              snumu := snu * umu.adj
              f[nu][ir] += cr * u[mu][ir] * snumu
              f[mu][ir] += cr * u[nu][ir] * snumu.adj
              ru[mu,nu][ir] += u[nu][ir].adj * u[mu][ir] * snu
        threadBarrier()
        sb[mu][nu].startSB(ru[mu,nu][ix])
        sb[nu][mu].startSB(ru[nu,mu][ix])
      toc("staple boundary")
      for ir in u[0]:
        for plane in planes:
          let (mu, nu) = plane
          if isLocal(sb[mu][nu],ir):
            var b: type(load1(u[0][0]))
            localSB(sb[mu][nu], ir, assign(b,it), ru[mu,nu][ix])
            f[mu][ir] += cr * b
          if isLocal(sb[nu][mu],ir):
            var b: type(load1(u[0][0]))
            localSB(sb[nu][mu], ir, assign(b,it), ru[nu,mu][ix])
            f[nu][ir] += cr * b
      toc("back rect local")
      for plane in planes:
        let (mu, nu) = plane
        var needBoundary = false
        boundaryWaitSB(sb[mu][nu]): needBoundary = true
        boundaryWaitSB(sb[nu][mu]): needBoundary = true
        if needBoundary:
          boundarySyncSB()
          for ir in lo:
            if not isLocal(sb[mu][nu],ir):
              var b: type(load1(u[0][0]))
              getSB(sb[mu][nu], ir, assign(b,it), ru[mu,nu][ix])
              f[mu][ir] += cr * b
            if not isLocal(sb[nu][mu],ir):
              var b: type(load1(u[0][0]))
              getSB(sb[nu][mu], ir, assign(b,it), ru[nu,mu][ix])
              f[nu][ir] += cr * b

proc fillStaplePlanes[V:static[int],T](uu: openArray[Field[V,T]], cs, ss: auto;
    stf, stu: StapleBatchFields[V,T]; planes: auto) =
  ## Fill one batch and start its backward exchanges. The caller completes
  ## those exchanges before rebinding or overwriting the scratch fields.
  mixin mul, load1, adj
  let u = cast[ptr cArray[Field[V,T]]](unsafeAddr(uu[0]))
  let lo = u[0].l
  threads:
    var umu, unu, umunu, unumu {.noinit.}: evalType(u[0][0])
    template fill(mu, nu, ir: untyped) =
      mul(umunu, umu, unu.adj)
      mul(stf[mu,nu][ir], u[nu][ir], umunu)
      mul(stf[nu,mu][ir], u[mu][ir], umunu.adj)
      mul(unumu, u[nu][ir].adj, u[mu][ir])
      mul(stu[mu,nu][ir], unumu, unu)
      mul(stu[nu,mu][ir], unumu.adj, umu)
    for ir in lo:
      for plane in planes:
        let (mu, nu) = plane
        if isLocal(cs[mu][nu],ir) and isLocal(cs[nu][mu],ir):
          localSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
          localSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
          fill(mu,nu,ir)
    for plane in planes:
      let (mu, nu) = plane
      var needBoundary = false
      boundaryWaitSB(cs[mu][nu]): needBoundary = true
      boundaryWaitSB(cs[nu][mu]): needBoundary = true
      if needBoundary:
        boundarySyncSB()
        for ir in lo:
          if not isLocal(cs[mu][nu],ir) or not isLocal(cs[nu][mu],ir):
            getSB(cs[mu][nu], ir, assign(umu,it), u[mu][ix])
            getSB(cs[nu][mu], ir, assign(unu,it), u[nu][ix])
            fill(mu,nu,ir)
      threadBarrier()
      ss[mu][nu].startSB(stu[mu,nu][ix])
      ss[nu][mu].startSB(stu[nu,mu][ix])

proc addStapleActions*[V:static[int],T](uu: openArray[Field[V,T]]; stf, stu: auto;
    ss: seq[seq[ShiftB[T]]]; planes: auto; rect: bool; sums: var seq[array[2,float]]) =
  ## sums[thread] += (sum Re tr U_mu stf[mu,nu], sum Re tr S stf for the
  ## backward staples S of stu shifted by ss when rect) over the thread's sites
  ## and the planes; ss is started. The caller sums over threads and ranks.
  mixin redot, load1
  let u = cast[ptr cArray[Field[V,T]]](unsafeAddr(uu[0]))
  let part = cast[ptr cArray[array[2,float]]](addr sums[0])
  let lo = u[0].l
  threads:
    var plaq, rec: float
    for ir in lo:
      for plane in planes:
        let (mu,nu) = plane
        plaq += simdSum(redot(u[mu][ir], stf[mu,nu][ir]))
        if rect:
          localSB(ss[mu][nu], ir,
            rec += simdSum(redot(it, stf[mu,nu][ir])), stu[mu,nu][ix])
          localSB(ss[nu][mu], ir,
            rec += simdSum(redot(it, stf[nu,mu][ir])), stu[nu,mu][ix])
    if rect:
      for plane in planes:
        let (mu,nu) = plane
        var needBoundary = false
        boundaryWaitSB(ss[mu][nu]): needBoundary = true
        boundaryWaitSB(ss[nu][mu]): needBoundary = true
        if needBoundary:
          boundarySyncSB()
          for ir in lo:
            if not isLocal(ss[mu][nu],ir):
              var b: type(load1(u[0][0]))
              getSB(ss[mu][nu], ir, assign(b,it), stu[mu,nu][ix])
              rec += simdSum(redot(b, stf[mu,nu][ir]))
            if not isLocal(ss[nu][mu],ir):
              var b: type(load1(u[0][0]))
              getSB(ss[nu][mu], ir, assign(b,it), stu[nu,mu][ix])
              rec += simdSum(redot(b, stf[nu,mu][ir]))
    part[threadNum][0] += plaq
    part[threadNum][1] += rec

proc stapleBatches*[V:static[int],T](cp, cr: float,
    uu: openArray[Field[V,T]], output: auto; derivative: static bool;
    accumulate = false): float =
  ## The action -(cp sum plaq + cr sum rect)/Nc, or its derivative into output,
  ## from staples of bounded batches of direction pairs; cp, cr are the raw
  ## plaquette and rectangle weights. Scratch is private to this call and
  ## released before returning; without a rect term the derivative needs no
  ## rectangle scratch.
  let u = cast[ptr cArray[Field[V,T]]](unsafeAddr(uu[0]))
  let lo = u[0].l
  let nd = lo.nDim
  let nc = float(u[0][0].ncols)
  if nd < 2:
    when derivative:
      if not accumulate:
        let f = cast[ptr cArray[type(output[0])]](unsafeAddr(output[0]))
        threads:
          for mu in 0..<nd: f[mu] := 0
    return 0.0
  let slots = min(temporaryStaplePairs, nd*(nd-1) div 2)
  let perPair = if derivative and cr != 0: 6 else: 4
  let scratch = newFieldArray(lo, type(uu[0]), perPair*slots)
  defer: release(scratch)
  var
    stf = StapleBatchFields[V,T](nd: nd, fields: newSeq[Field[V,T]](nd*nd))
    stu = StapleBatchFields[V,T](nd: nd, fields: newSeq[Field[V,T]](nd*nd))
  let cs = startCornerShifts(uu)
  let ss = newLinkShifts(lo, T, -1)
  when derivative:
    var ru = StapleBatchFields[V,T](nd: nd, fields: newSeq[Field[V,T]](nd*nd))
    let sf = if cr != 0: newLinkShifts(lo, T, 1) else: newSeq[seq[ShiftB[T]]]()
    let sb = if cr != 0: newLinkShifts(lo, T, -1) else: newSeq[seq[ShiftB[T]]]()
  else:
    var sums = newSeq[array[2,float]](getMaxThreads())
  let planes = staplePlanes(nd)
  for first in countup(0, planes.high, slots):
    let batch = planes[first..<min(first+slots,planes.len)]
    for slot, plane in batch:
      let (mu,nu) = plane
      stf.fields[mu*nd+nu] = scratch[2*slot]
      stf.fields[nu*nd+mu] = scratch[2*slot+1]
      stu.fields[mu*nd+nu] = scratch[2*slots+2*slot]
      stu.fields[nu*nd+mu] = scratch[2*slots+2*slot+1]
      when derivative:
        if cr != 0:
          ru.fields[mu*nd+nu] = scratch[4*slots+2*slot]
          ru.fields[nu*nd+mu] = scratch[4*slots+2*slot+1]
    fillStaplePlanes(uu, cs, ss, stf, stu, batch)
    when derivative:
      addStaplePlanes(uu, output, cp/nc, cr/nc,
        stf, stu, ru, cs, ss, sf, sb, batch, overwrite=not accumulate)
    else:
      addStapleActions(uu, stf, stu, ss, batch, cr != 0, sums)
  when not derivative:
    var a = [0.0, 0.0]
    for part in sums:
      a[0] += part[0]
      a[1] += part[1]
    rankSum(a)
    result = -(cp*a[0]+cr*a[1])/nc

when isMainModule:
  import qex
  import physics/qcdTypes

  qexInit()
  #var defaultGaugeFile = "l88.scidac"
  #let defaultLat = @[2,2,2,2]
  let defaultLat = @[8,8,8,8]
  defaultSetup()
  g.unit
  #g.random

  proc test(g: auto) =
    tic("test")
    var cs = startCornerShifts(g)
    toc("startCornerShifts")
    var (stf,stu,ss) = makeStaples(g, cs)
    toc("makeStaples")

  test(g)
  resetTimers()
  test(g)
  echoTimers()

  qexFinalize()
