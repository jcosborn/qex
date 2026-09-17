import qex, comms/halo
import std/math
export halo

const
  sqrt_1_2* = 0.707106781186547524400844362104849039    # sqrt(1/2)
  sqrt_1_2pi* = 0.398942280401432677939946059934381868  # sqrt(1/(2*pi))

proc realField*[T: SomeFloat](lo: Layout[VLEN]; precision: typedesc[T]): auto =
  when T is float32: lo.RealS()
  else: lo.RealD()

type RealField*[T: SomeFloat] = typeof(realField(default(Layout[VLEN]),T))
type SiteMask* = RealField[float32]
  ## Existing lattice field storage; zero excludes a site, nonzero selects it.

proc selected*(mask: SiteMask; site: int): bool {.inline.} =
  if mask == nil: return true
  var x: float32
  x := mask{site}
  x != 0

template store(f, e, mask, R, value: untyped) =
  ## f[e] := value on the lanes mask selects; nil selects every lane. A SIMD
  ## site whose lanes are all selected or all excluded is stored whole.
  if mask == nil:
    f[e] := value
  else:
    let m = abs(mask[e])
    if simdMin(m) != 0:
      f[e] := value
    elif simdMax(m) != 0:
      let r = value
      for k in 0..<simdLength(m):
        var w: float32
        w := m[asSimd(k)]
        if w != 0:
          var v: R
          v := r[asSimd(k)]
          f{e*simdLength(m)+k} := v

type
  ConvParams*[T] = object
    cin*, cout*: int
    kernel*: seq[int]
    offsets*: seq[seq[int32]]
    weights*, bias*: seq[T]
  ConvWorkspace*[L,F,E] = ref object
    layout*: HaloLayout[L]
    map*: HaloMap[L]
    halo*: seq[Halo[L,F,E]]
    offsets*: seq[seq[int32]]
    index*: seq[int32]

proc validate[T](p: ConvParams[T], nd: int, weights = true) =
  if p.cin <= 0 or p.cout <= 0 or p.offsets.len == 0:
    raise newException(ValueError, "convolution requires channels and taps")
  if weights and p.weights.len != p.cin*p.cout*p.offsets.len:
    raise newException(ValueError, "convolution weight length does not match channels and taps")
  if p.bias.len != 0 and p.bias.len != p.cout:
    raise newException(ValueError, "convolution bias length does not match output channels")
  for off in p.offsets:
    if off.len != nd:
      raise newException(ValueError, "convolution offset dimension does not match layout")

proc kernelOffsets*(shape: openArray[int]): seq[seq[int32]] =
  ## Taps of a dense kernel with odd extents, last spatial dimension fastest.
  if shape.len == 0:
    raise newException(ValueError, "convolution requires a spatial kernel")
  var n = 1
  for k in shape:
    if k <= 0 or k mod 2 != 1:
      raise newException(ValueError, "convolution kernel extents must be positive and odd")
    n *= k
  result = newSeq[seq[int32]](n)
  for t in 0..<n:
    var j = t
    result[t] = newSeq[int32](shape.len)
    for d in countdown(shape.len-1, 0):
      result[t][d] = int32(j mod shape[d] - shape[d] div 2)
      j = j div shape[d]

proc convParams*[T](cin, cout: int; shape: openArray[int]; weights: seq[T]; bias: seq[T] = @[]): ConvParams[T] =
  ## Weights are [output,input,tap], with the last spatial dimension fastest.
  result.cin = cin
  result.cout = cout
  result.kernel = @shape
  result.weights = weights
  result.bias = bias
  result.offsets = kernelOffsets(shape)
  result.validate(shape.len)

proc convWorkspace*[F,T](proto: F; p: ConvParams[T]): auto =
  ## Halo buffers, exchange map and tap table for inputs shaped like proto.
  ## Every convolution rebinds the halo fields, so one workspace serves any
  ## inputs of that layout and channel count.
  let lo = proto.l
  p.validate(lo.nDim, false)
  let hl = lo.haloLayout(p.offsets)
  type E = eval(F.type.index(int))
  type L = type(lo)
  var ws: ConvWorkspace[L,F,E]
  ws.new
  ws.layout = hl
  ws.offsets = p.offsets
  ws.halo = newSeq[Halo[L,F,E]](p.cin)
  if hl.nExt > hl.nOut:
    ws.map = hl.haloMap(lo.comm, p.offsets)
  for i in 0..<p.cin:
    ws.halo[i] = hl.makeHalo(proto)
  ws.index = hl.haloIndex(p.offsets)
  ws

proc validateFields[F](dst, src: openArray[F]; mask: SiteMask) =
  if dst.len == 0 or src.len == 0:
    raise newException(ValueError, "NN operations require channels")
  let lo = src[0].l
  for f in src:
    if f.l != lo: raise newException(ValueError, "NN source layouts differ")
  for f in dst:
    if f.l != lo: raise newException(ValueError, "NN destination layout differs from source")
  if mask != nil and mask.l != lo:
    raise newException(ValueError, "NN mask layout differs from source")

proc conv*[L,F,E,T](dst, src: seq[F]; p: ConvParams[T]; ws: ConvWorkspace[L,F,E]; sub = "all"; addBias = true) =
  ## Call outside a threads block. Workspace buffers are refreshed on every call.
  ## z_o(x) = sum_i,t W_o,i,t x_i(x+offset_t), followed by bias. A masked
  ## output composes with maskedCopy.
  validateFields(dst, src, nil)
  if dst.len != p.cout or src.len != p.cin:
    raise newException(ValueError, "convolution channel count does not match parameters")
  p.validate(src[0].l.nDim)
  if ws.layout.lo != src[0].l or ws.halo.len != p.cin or ws.offsets != p.offsets:
    raise newException(ValueError, "convolution workspace does not match parameters and layout")
  for d in dst:
    for s in src:
      if d == s:
        raise newException(ValueError, "convolution destination must not alias an input")
  for i in 0..<p.cin:
    ws.halo[i].field = src[i]
    if ws.layout.nExt > ws.layout.nOut:
      ws.halo[i].update(ws.map, ws.layout.lo.comm)
  let nt = p.offsets.len
  threads:
    for o in 0..<p.cout:
      for x in dst[o][sub]:
        var z: E
        z := 0
        for i in 0..<p.cin:
          for t in 0..<nt:
            let j = ws.index[x*nt+t]
            z += p.weights[(o*p.cin+i)*nt+t] * ws.halo[i][j]
        if addBias and p.bias.len != 0:
          z += p.bias[o]
        dst[o][x] := z

template mapFields(dst, src, sub, mask, expression: untyped) =
  ## expression of the SIMD site value x, stored on the selected lanes of sub.
  validateFields(dst, src, mask)
  if dst.len != src.len:
    raise newException(ValueError, "pointwise channel counts differ")
  type R {.inject.} = numberType(typeof(src[0]))
  threads:
    for c in 0..<src.len:
      for e in dst[c][sub]:
        var x {.inject.}: evalType(src[c][e])
        x := src[c][e]
        store(dst[c], e, mask, R, expression)

proc gelu*[F](dst, src: seq[F]; sub = "all"; mask: SiteMask = nil) =
  ## Exact GELU: x Phi(x). Supports elementwise in-place evaluation.
  mapFields(dst, src, sub, mask): (R(0.5)*x) * erfc(-x*R(sqrt_1_2))

proc arctan*[F](dst, src: seq[F]; sub = "all"; mask: SiteMask = nil) =
  mapFields(dst, src, sub, mask): arctan(x)

proc divide*[F; T: SomeFloat](dst, src: seq[F]; divisor: T; sub = "all"; mask: SiteMask = nil) =
  ## Round the divisor to the field precision; each call performs one division.
  mapFields(dst, src, sub, mask): x / R(divisor)

proc exp*[F](dst, src: seq[F]; sub = "all"; mask: SiteMask = nil) =
  mapFields(dst, src, sub, mask): exp(x)

proc ln*[F](dst, src: seq[F]; sub = "all"; mask: SiteMask = nil) =
  mapFields(dst, src, sub, mask): ln(x)

proc erfc*[F](dst, src: seq[F]; sub = "all"; mask: SiteMask = nil) =
  mapFields(dst, src, sub, mask): erfc(x)

proc clipMin*[F; T: SomeFloat](dst, src: seq[F]; floor: T; sub = "all"; mask: SiteMask = nil) =
  ## max(x, floor); clipSlope is its derivative in the JAX maximum convention.
  var f: eval(F.type.index(int))
  f := numberType(F)(floor)
  mapFields(dst, src, sub, mask): max(x, f)

proc slope[E; R: SomeFloat](x: E; f: R): E =
  ## 0 below f, 1/2 at equality, 1 above, lane by lane.
  when simdLength(E) == 1:
    result = (if x < f: R(0) elif x == f: R(0.5) else: R(1))
  else:
    for k in 0..<simdLength(E):
      var v: R
      v := x[asSimd(k)]
      result[asSimd(k)] = (if v < f: R(0) elif v == f: R(0.5) else: R(1))

proc clipSlope*[F; T: SomeFloat](dst, src: seq[F]; floor: T; sub = "all"; mask: SiteMask = nil) =
  ## 0 below floor, 1/2 at equality, 1 above.
  let f = numberType(F)(floor)
  mapFields(dst, src, sub, mask): slope(x, f)

proc divide*[F](dst, x, y: seq[F]; sub = "all"; mask: SiteMask = nil) =
  ## Elementwise x/y in the field precision.
  validateFields(dst, x, mask)
  validateFields(dst, y, nil)
  if dst.len != x.len or dst.len != y.len:
    raise newException(ValueError, "division channel counts differ")
  type R = numberType(F)
  threads:
    for c in 0..<dst.len:
      for e in dst[c][sub]:
        store(dst[c], e, mask, R, x[c][e]/y[c][e])

proc redot*[F](x, y: seq[F]): float =
  ## sum_c sum_s x_c(s) y_c(s) over the global lattice, accumulated in double.
  validateFields(x, y, nil)
  if x.len != y.len:
    raise newException(ValueError, "dot product channel counts differ")
  var total = 0.0
  threads:
    var acc: DLatticeRealV.T
    acc := 0
    for c in 0..<x.len:
      for e in x[c]:
        var a, b: DLatticeRealV.T
        a := x[c][e]
        b := y[c][e]
        acc += a*b
    var s = simdSum(acc)
    x[0].l.threadRankSum(s)
    threadSingle: total = s
  total

proc scale*[F,T](dst, src: seq[F]; scales: seq[T]; sub = "all"; mask: SiteMask = nil) =
  validateFields(dst, src, mask)
  if dst.len != src.len or scales.len != src.len:
    raise newException(ValueError, "channel scale count differs from field channels")
  type R = numberType(F)
  threads:
    for c in 0..<src.len:
      for e in dst[c][sub]:
        store(dst[c], e, mask, R, src[c][e] * scales[c])

proc bias*[F,T](dst, src: seq[F]; biases: seq[T]; sub = "all"; mask: SiteMask = nil) =
  validateFields(dst, src, mask)
  if dst.len != src.len or biases.len != src.len:
    raise newException(ValueError, "channel bias count differs from field channels")
  type R = numberType(F)
  threads:
    for c in 0..<src.len:
      for e in dst[c][sub]:
        store(dst[c], e, mask, R, src[c][e] + biases[c])

proc maskedCopy*[F](dst, src: seq[F]; mask: SiteMask; sub = "all") =
  validateFields(dst, src, mask)
  if dst.len != src.len:
    raise newException(ValueError, "masked copy channel counts differ")
  type R = numberType(F)
  threads:
    for c in 0..<src.len:
      for e in dst[c][sub]:
        store(dst[c], e, mask, R, src[c][e])

proc convVjp*[L,F,E,T](dx, dy: seq[F]; p: ConvParams[T]; ws: ConvWorkspace[L,F,E]; sub = "all") =
  ## dx_i(x+tap) += W_o,i,tap dy_o(x), followed by reverse halo accumulation.
  ## Each thread owns complete input channels, including their halo gradients.
  ## Writes dx from zero; workspace halos are rebound and reusable by conv.
  ## A masked output's seed is maskVjp of dy.
  validateFields(dx, dy, nil)
  if dx.len != p.cin or dy.len != p.cout:
    raise newException(ValueError, "convolution VJP channel count does not match parameters")
  let lo = dy[0].l
  p.validate(lo.nDim)
  if ws.layout.lo != lo or ws.halo.len != p.cin or ws.offsets != p.offsets:
    raise newException(ValueError, "convolution VJP workspace does not match parameters and layout")
  for a in dx:
    for b in dy:
      if a == b:
        raise newException(ValueError, "convolution VJP destination must not alias its seed")
  for i in 0..<p.cin: ws.halo[i].field = dx[i]
  let sel = lo.getSubset(sub)
  let nt = p.offsets.len
  threads:
    var i = threadNum
    while i < p.cin:
      let h = ws.halo[i]
      for x in 0..<h.nExt: h[x] := 0
      for o in 0..<p.cout:
        for x in sel.lowOuter..<sel.highOuter:
          var v: E
          v := dy[o][x]
          for t in 0..<nt:
            let j = ws.index[x*nt+t]
            h[j] += p.weights[(o*p.cin+i)*nt+t] * v
      i += numThreads
  if ws.layout.nExt > ws.layout.nOut:
    for i in 0..<p.cin: ws.halo[i].updateRev(ws.map, lo.comm)

proc convWeightVjp*[L,F,E,T](dw: var seq[T]; x, dy: seq[F]; p: ConvParams[T]; ws: ConvWorkspace[L,F,E]; sub = "all") =
  ## dw_o,i,t = sum_s dy_o(s) x_i(s+offset_t), globally replicated.
  ## Threads own parameter entries; the rank reduction runs once after joining.
  validateFields(x,dy,nil)
  let lo = x[0].l
  p.validate(lo.nDim, false)
  let nt = p.offsets.len
  if x.len != p.cin or dy.len != p.cout or dw.len != p.cout*p.cin*nt:
    raise newException(ValueError, "convolution weight VJP shape differs from parameters")
  if ws.layout.lo != lo or ws.halo.len != p.cin or ws.offsets != p.offsets:
    raise newException(ValueError, "convolution weight VJP workspace does not match parameters and layout")
  for i in 0..<p.cin:
    ws.halo[i].field = x[i]
    if ws.layout.nExt > ws.layout.nOut:
      ws.halo[i].update(ws.map, lo.comm)
  let sel = lo.getSubset(sub)
  let dst = cast[ptr UncheckedArray[T]](dw[0].addr)
  let n = dw.len
  threads:
    var k = threadNum
    while k < n:
      let t = k mod nt
      let i = (k div nt) mod p.cin
      let o = k div (p.cin*nt)
      var acc: E
      acc := 0
      for s in sel.lowOuter..<sel.highOuter:
        var v: E
        v := dy[o][s]
        acc += v*ws.halo[i][ws.index[s*nt+t]]
      dst[k] = T(simdSum(acc))
      k += numThreads
  lo.comm.rankSum(dw)

proc channelSum*[F,T](dst: var seq[T]; src: seq[F]) =
  ## dst_c = sum_s src_c(s), with one global reduction of the channel array.
  validateFields(src,src,nil)
  if dst.len != src.len:
    raise newException(ValueError, "channel sum count differs from fields")
  type E = eval(F.type.index(int))
  let lo = src[0].l
  let data = cast[ptr UncheckedArray[T]](dst[0].addr)
  threads:
    var c = threadNum
    while c < src.len:
      var acc: E
      acc := 0
      for s in 0..<lo.nSitesOuter: acc += src[c][s]
      data[c] = T(simdSum(acc))
      c += numThreads
  lo.comm.rankSum(dst)

proc broadcast*[F,T](dst: seq[F]; values: seq[T]) =
  ## dst_c(s) = values_c; the parameter array is already replicated.
  validateFields(dst,dst,nil)
  if dst.len != values.len:
    raise newException(ValueError, "broadcast count differs from fields")
  threads:
    for c in 0..<dst.len: dst[c] := values[c]

template mapVjp(dx, dy, sub, mask, passthrough, body: untyped) =
  ## body maps the SIMD seed b on sub's sites that mask selects; elsewhere dx
  ## keeps the seed when passthrough and is zero otherwise.
  validateFields(dx, dy, mask)
  if dx.len != dy.len:
    raise newException(ValueError, "pointwise VJP channel counts differ")
  type R {.inject.} = numberType(typeof(dy[0]))
  let sel = dy[0].l.getSubset(sub)
  threads:
    for c {.inject.} in 0..<dx.len:
      for e {.inject.} in dx[c]:
        var b {.inject.}: evalType(dy[c][e])
        b := dy[c][e]
        if e >= sel.lowOuter and e < sel.highOuter:
          if mask == nil:
            body
            dx[c][e] := b
          else:
            let keep = b
            body
            if passthrough: dx[c][e] := keep
            else: dx[c][e] := 0
            store(dx[c], e, mask, R, b)
        else:
          if not passthrough: b := 0
          dx[c][e] := b

proc geluVjp*[F](dx, x, dy: seq[F]; sub = "all"; mask: SiteMask = nil; passthrough = false) =
  ## Masked fill has zero input derivative; passthrough preserves its complement.
  validateFields(dx, x, mask)
  if dx.len != x.len: raise newException(ValueError, "GELU VJP channel counts differ")
  mapVjp(dx, dy, sub, mask, passthrough):
    var a: typeof(b)
    a := x[c][e]
    let phi = R(0.5)*erfc(-a*R(sqrt_1_2))
    let pdf = R(sqrt_1_2pi)*exp((-R(0.5)*a)*a)
    b *= phi+a*pdf

proc arctanVjp*[F](dx, x, dy: seq[F]; sub = "all"; mask: SiteMask = nil; passthrough = false) =
  validateFields(dx, x, mask)
  if dx.len != x.len: raise newException(ValueError, "arctan VJP channel counts differ")
  mapVjp(dx, dy, sub, mask, passthrough):
    var a: typeof(b)
    a := x[c][e]
    b = b/(R(1)+a*a)

proc divideVjp*[F; T: SomeFloat](dx, dy: seq[F]; divisor: T; sub = "all"; mask: SiteMask = nil; passthrough = false) =
  mapVjp(dx, dy, sub, mask, passthrough): b = b/R(divisor)

proc scaleVjp*[F,T](dx, dy: seq[F]; scales: seq[T]; sub = "all"; mask: SiteMask = nil; passthrough = false) =
  if scales.len != dy.len:
    raise newException(ValueError, "scale VJP channel counts differ")
  mapVjp(dx, dy, sub, mask, passthrough): b *= R(scales[c])

proc maskVjp*[F](dx, dy: seq[F]; sub = "all"; mask: SiteMask = nil; complement = false) =
  ## Source/bias-input cotangent selects the written sites. complement=true
  ## gives the old destination/fill cotangent of a masked assignment.
  mapVjp(dx, dy, sub, mask, complement):
    if complement: b := 0
