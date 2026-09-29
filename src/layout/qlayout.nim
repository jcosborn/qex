import strutils
import base
import layoutTypes

proc `$`(a: ptr cArray): string =
  result.add $a[0]
  for i in 1..3:
    result.add(" " & $a[i])

proc layoutSetupQ*(l: var LayoutQ) =
  var nd = l.nDim
  l.outerGeom = cast[type(l.outerGeom)](alloc(nd * sizeof(cint)))
  l.localGeom = cast[type(l.localGeom)](alloc(nd * sizeof(cint)))
  var
    pvol = 1
    lvol = 1
    ovol = 1
    icb = 0
    icbd = -1
  for i in 0..<nd:
    l.localGeom[i] = l.physGeom[i] div l.rankGeom[i]
    l.outerGeom[i] = l.localGeom[i] div l.innerGeom[i]
    pvol = pvol * l.physGeom[i]
    lvol = lvol * l.localGeom[i]
    ovol = ovol * l.outerGeom[i]
    if l.innerGeom[i] > 1 and (l.outerGeom[i] and 1) == 1: inc(icb)
    if l.innerGeom[i] == 1 and (l.outerGeom[i] and 1) == 0: icbd = i
  if icb == 0:
    icbd = 0
  else:
    if icbd < 0:
      if l.myrank == 0:
        echo "not enough 2\'s in localGeom"
        echo "physGeom: ", l.physGeom
        echo "rankGeom: ", l.rankGeom
        echo "localGeom: ", l.localGeom
        echo "outerGeom: ", l.outerGeom
        echo "innerGeom: ", l.innerGeom
      quit(-1)
    icb = l.outerGeom[icbd] div 2
    if (icb and 1) == 0:
      if l.myrank == 0:
        echo "error in cb choice"
        echo "physGeom: ", l.physGeom
        echo "rankGeom: ", l.rankGeom
        echo "localGeom: ", l.localGeom
        echo "outerGeom: ", l.outerGeom
        echo "innerGeom: ", l.innerGeom
        echo "innerCb: ", icb
        echo "innerCbDir: ", icbd
      quit(-1)
  l.physVol = pvol
  l.nSites = lvol
  l.nOdd = lvol div 2
  l.nEven = lvol - l.nOdd
  l.nSitesOuter = ovol
  l.nOddOuter = ovol div 2
  l.nEvenOuter = ovol - l.nOddOuter
  l.nSitesInner = int32(l.nSites div l.nSitesOuter)
  l.innerCb = int32 icb
  l.innerCbDir = int32 icbd
  if l.myrank == 0:
    echo "#innerCb: ", icb
    echo "#innerCbDir: ", icbd


proc lex_x*(x: var openArray; ll: SomeInteger; s: openArray; ndim: SomeInteger) =
  var l = ll
  for i in 0..<ndim:
    x[i] = l mod s[i]
    l = l div s[i]

proc lexr_x*[T](x: var openArray[T]; ll: SomeInteger;
             s: ptr cArray[cint]; ndim: cint) =
  var l = ll
  for i in countdown(ndim-1, 0):
    x[i] = l mod s[i]
    l = l div s[i]

# x[0] is fastest
proc lex_i*[X,S:UncheckedArray[SomeInteger],N:SomeInteger](
  x: ptr X, s: ptr S, d: ptr UncheckedArray[int32]; ndim: N): int =
  var l = 0
  #var i: cint = ndim - 1
  #while i >= 0:
  for i in countdown(ndim-1,0):
    var xx = x[i]
    if not d.isNil: xx = xx div d[i]
    l = l * s[i] + (xx mod s[i])
    #dec(i)
  return l

# x[0] is slowest
proc lexr_i*[X,S,D:UncheckedArray[SomeInteger],N:SomeInteger](
  x: ptr X, s: ptr S, d: ptr D; ndim: N): int =
  var l = 0
  #var i = 0
  #while i < ndim:
  for i in 0..<ndim:
    var xx = x[i]
    if not d.isNil: xx = xx div d[i]
    l = l * s[i] + (xx mod s[i])
    #inc(i)
  return l

#proc layoutLocalIndexQ*[T](l: LayoutQ; coords: var openArray[T]): int32 =

template layoutIndexImpl(l: LayoutQ; li: var LayoutIndexQ; crd: untyped) =
  # x = crd(d), the global coordinate in direction d:
  #   rank  (x div localGeom) mod rankGeom -> ri, direction 0 slowest
  #   lane  k = (x div outerGeom) mod innerGeom -> ii, direction 0 fastest
  #   outer o = x mod outerGeom -> oi, direction 0 fastest
  # ib = sum k*outerGeom; odd ib shifts o[innerCbDir] by innerCb mod outerGeom.
  # index = ((oi + (sum x odd)*nSitesOuter) div 2)*nSitesInner + ii
  var ri, ii, oi, ib, p, ocb, mcb = 0
  var mi, mo = 1
  for d in 0..<l.nDim.int:
    let x = int(crd(d))
    let og = int(l.outerGeom[d])
    let o = x mod og
    let k = (x div og) mod l.innerGeom[d]
    ri = ri*l.rankGeom[d] + (x div l.localGeom[d]) mod l.rankGeom[d]
    ii += k*mi
    mi *= l.innerGeom[d]
    if d == l.innerCbDir:
      ocb = o
      mcb = mo
    oi += o*mo
    mo *= og
    ib += k*og
    p += x
  if (ib and 1) != 0:
    oi += ((ocb + l.innerCb) mod l.outerGeom[l.innerCbDir] - ocb)*mcb
  if (p and 1) != 0: oi += l.nSitesOuter
  li.rank = int32 ri
  li.index = int32((oi div 2)*l.nSitesInner + ii)

proc layoutIndexQ*[T](l: LayoutQ; li: var LayoutIndexQ; coords: openArray[T]) =
  template crd(d: int): untyped = coords[d]
  layoutIndexImpl(l, li, crd)

proc layoutIndexQ*(l: LayoutQ; li: var LayoutIndexQ; coords: seq[seq[int16]]; i: int) =
  ## Coordinates coords[d][i], as in the Layout.coords table.
  template crd(d: int): untyped = coords[d][i]
  layoutIndexImpl(l, li, crd)

proc layoutCoordQ*[T](l: ptr LayoutQ; coords: var openArray[T];
                      li: ptr LayoutIndexQ) =
  # coords[i] = localGeom[i]*rank[i] + outerGeom[i]*lane[i] + outer[i]
  # localGeom[i] = innerGeom[i]*outerGeom[i], so the first two terms are
  # multiples of outerGeom[i] and outer[i] = coords[i] mod outerGeom[i].
  var nd = l.nDim
  var r = li.rank
  for i in countdown(nd-1, 0):
    coords[i] = l.localGeom[i] * (r mod l.rankGeom[i])
    r = r div l.rankGeom[i]
  var p = 0
  var ll = li.index mod l.nSitesInner
  var ib = 0
  for i in 0..<nd:
    var w = l.innerGeom[i]
    var wl = l.outerGeom[i]
    var k = ll mod w
    coords[i] += k * wl
    inc(p, coords[i])
    ll = ll div w
    inc(ib, k * wl)
  ib = ib and 1
  var ii = li.index div l.nSitesInner
  if ii >= l.nEvenOuter:
    dec(ii, l.nEvenOuter)
    inc(p)
  ii = ii * 2
  for i in 0..<nd:
    var wl = l.outerGeom[i]
    var k = ii mod wl
    if i == l.innerCbDir: k = (k + l.innerCb * ib).int32 mod wl
    coords[i] += k
    inc(p, k)
    ii = ii div wl
  if (p and 1) != 0:
    for i in 0..<nd:
      var wl: cint = l.outerGeom[i]
      var k = coords[i] mod wl
      let b = coords[i] - k
      if i == l.innerCbDir: k = int32(k + l.innerCb * ib) mod wl
      inc(k)
      if k >= wl:
        k = 0
        if i == l.innerCbDir: k = int32(k + l.innerCb * ib) mod wl
        coords[i] = b + k
      else:
        if i == l.innerCbDir: k = int32(k + l.innerCb * ib) mod wl
        coords[i] = b + k
        break
  var li2: LayoutIndexQ
  layoutIndexQ(l[], li2, coords)
  if li.rank != li2.rank or li.index != li2.index:
    echo "error: bad coord:"
    echo " $#,$# -> $# $# $# $# -> $#,$#"%[$li.rank,$li.index,$coords[0],
           $coords[1], $coords[2], $coords[3], $li2.rank, $li2.index]
    quit(-1)

#[
proc layoutShift*(l: LayoutQ; li: LayoutIndexQ; li2: LayoutIndexQ; 
                  disp: openarray) = 
  var nd: cint = l.nDim
  var x: array[nd, cint]
  layoutCoord(l, x, li2)
  var i: cint = 0
  while i < nd: 
    x[i] = (x[i] + disp[i] + l.physGeom[i]) mod l.physGeom[i]
    inc(i)
  layoutIndex(l, li, x)
]#
