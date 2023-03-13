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

template `&`[T](x: openArray[T]): untyped = cast[ptr cArray[T]](unsafeaddr x[0])

#proc layoutLocalIndexQ*[T](l: LayoutQ; coords: var openArray[T]): int32 =

proc layoutIndexQ*[T](l: LayoutQ; li: var LayoutIndexQ;
                      coords: var openArray[T]) =
  var nd = l.nDim
  var ri = lexr_i(&coords, l.rankGeom, l.localGeom, nd)
  #var ri = lex_i(coords, l.rankGeom, l.localGeom, nd)
  var ii = lex_i(&coords, l.innerGeom, l.outerGeom, nd)
  var ib = 0
  for i in 0..<nd:
    var xi = coords[i] div l.outerGeom[i]
    var xli = xi mod l.innerGeom[i]
    inc(ib, xli * l.outerGeom[i])
  ib = ib and 1
  coords[l.innerCbDir] += int32(l.innerCb * ib)
  var oi = lex_i(&coords, l.outerGeom, nil, nd)
  coords[l.innerCbDir] -= int32(l.innerCb * ib)
  var p = 0
  for i in 0..<nd:
    inc(p, coords[i])
  var oi2 = oi div 2
  if (p and 1) != 0: oi2 = (oi + l.nSitesOuter).int32 div 2
  li.rank = int32 ri
  li.index = int32 oi2 * l.nSitesInner + ii

proc layoutCoordQ*[T](l: ptr LayoutQ; coords: var openArray[T];
                      li: ptr LayoutIndexQ) =
  var nd = l.nDim
  var cr = newSeq[cint](nd)
  lexr_x(cr, li.rank, l.rankGeom, nd)
  #lex_x(cr, li.rank, l.rankGeom, nd)
  var p = 0
  var ll = li.index mod l.nSitesInner
  var ib = 0
  for i in 0..<nd:
    var w = l.innerGeom[i]
    var wl = l.outerGeom[i]
    var k = ll mod w
    var c = l.localGeom[i] * cr[i] + k * wl
    cr[i] = c
    #printf("cr[%i]: %i\n", i, c);
    inc(p, c)
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
    coords[i] = k
    #printf("coords[%i]: %i\n", i, k);
    inc(p, k)
    ii = ii div wl
  if (p and 1) != 0:
    for i in 0..<nd:
      var wl: cint = l.outerGeom[i]
      if i == l.innerCbDir: coords[i] = int32(coords[i] + l.innerCb * ib) mod wl
      inc(coords[i])
      if coords[i] >= wl:
        coords[i] = 0
        if i == l.innerCbDir: coords[i] = int32(coords[i] + l.innerCb * ib) mod wl
      else:
        if i == l.innerCbDir: coords[i] = int32(coords[i] + l.innerCb * ib) mod wl
        break
  for i in 0..<nd:
    coords[i] += cr[i]
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
