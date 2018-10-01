import strutils
import base
import layoutTypes
#const 
#  myalloc* = malloc
#template PRINTV*(s, f, v, n: expr): stmt = 
#  while true: 
#  printf(s)
#  var _i: cint = 0
#  while _i < n: 
#    printf(" ", f, (v)[_i])
#    inc(_i)
#  printf("\x0A")
#  if not 0: break 

proc `$`(a: ptr cArray[system.cint]): string =
  result.add $a[0]
  for i in 1..3:
    result.add(" " & $a[i])

proc layoutSetupQ*(l: ptr LayoutQ) =
  var nd: cint = l.nDim
  l.outerGeom = cast[type(l.outerGeom)](alloc(nd * sizeof(cint)))
  l.localGeom = cast[type(l.localGeom)](alloc(nd * sizeof(cint)))
  var
    pvol: cint = 1
    lvol: cint = 1
    ovol: cint = 1
    icb: cint = 0
    icbd: cint = - 1
  var i: cint = 0
  while i < nd:
    l.localGeom[i] = l.physGeom[i] div l.rankGeom[i]
    l.outerGeom[i] = l.localGeom[i] div l.innerGeom[i]
    pvol = pvol * l.physGeom[i]
    lvol = lvol * l.localGeom[i]
    ovol = ovol * l.outerGeom[i]
    if l.innerGeom[i] > 1 and (l.outerGeom[i] and 1) == 1: inc(icb)
    if l.innerGeom[i] == 1 and (l.outerGeom[i] and 1) == 0: icbd = i
    inc(i)
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
  l.nSitesInner = l.nSites div l.nSitesOuter
  l.innerCb = icb
  l.innerCbDir = icbd
  if l.myrank == 0:
    echo "#innerCb: ", icb
    echo "#innerCbDir: ", icbd

proc lex_x*(x: var openArray[cint]; ll: cint;
            s: ptr cArray[cint]; ndim: cint) =
  var i: cint = 0
  var l = ll
  while i < ndim:
    x[i] = l mod s[i]
    l = l div s[i]
    inc(i)

# x[0] is fastest
proc lex_i*(x: ptr cArray[cint]; s: ptr cArray[cint];
            d: ptr cArray[cint]; ndim: cint): cint =
  var l: cint = 0
  var i: cint = ndim - 1
  while i >= 0:
    var xx: cint = x[i]
    if not d.isNil: xx = xx div d[i]
    l = l * s[i] + (xx mod s[i])
    dec(i)
  return l

when false:
  # x[0] is slowest
  proc lexr_i*(x: ptr cint; s: ptr cint; d: ptr cint; ndim: cint): cint = 
    var l: cint = 0
    var i: cint = 0
    while i < ndim: 
      var xx: cint = x[i]
      if d: xx = xx div d[i]
      l = l * s[i] + (xx mod s[i])
      inc(i)
    return l

proc layoutIndexQ*(l: ptr LayoutQ; li: ptr LayoutIndexQ;
                   coords: ptr cArray[cint]) =
  var nd: cint = l.nDim
  var ri: cint = lex_i(coords, l.rankGeom, l.localGeom, nd)
  var ii: cint = lex_i(coords, l.innerGeom, l.outerGeom, nd)
  var ib: cint = 0
  var i: cint = 0
  while i < nd:
    var xi: cint = coords[i] div l.outerGeom[i]
    var li: cint = xi mod l.innerGeom[i]
    inc(ib, li * l.outerGeom[i])
    inc(i)
  ib = ib and 1
  inc(coords[l.innerCbDir], l.innerCb * ib)
  var oi: cint = lex_i(coords, l.outerGeom, nil, nd)
  dec(coords[l.innerCbDir], l.innerCb * ib)
  var p: cint = 0
  i=0
  while i < nd:
    inc(p, coords[i])
    inc(i)
  var oi2: cint = oi div 2
  if (p and 1) != 0: oi2 = (oi + l.nSitesOuter) div 2
  li.rank = ri
  li.index = oi2 * l.nSitesInner + ii

proc layoutCoordQ*(l: ptr LayoutQ; coords: ptr cArray[cint];
                   li: ptr LayoutIndexQ) =
  var nd: cint = l.nDim
  var cr = newSeq[cint](nd)
  lex_x(cr, li.rank, l.rankGeom, nd)
  var p: cint = 0
  var ll: cint = li.index mod l.nSitesInner
  var ib: cint = 0
  var i: cint = 0
  while i < nd:
    var w: cint = l.innerGeom[i]
    var wl: cint = l.outerGeom[i]
    var k: cint = ll mod w
    var c: cint = l.localGeom[i] * cr[i] + k * wl
    cr[i] = c
    #printf("cr[%i]: %i\n", i, c);
    inc(p, c)
    ll = ll div w
    inc(ib, k * wl)
    inc(i)
  ib = ib and 1
  var ii: cint = li.index div l.nSitesInner
  if ii >= l.nEvenOuter:
    dec(ii, l.nEvenOuter)
    inc(p)
  ii = ii * 2
  i=0
  while i < nd:
    var wl: cint = l.outerGeom[i]
    var k: cint = ii mod wl
    if i == l.innerCbDir: k = (k + l.innerCb * ib) mod wl
    coords[i] = k
    #printf("coords[%i]: %i\n", i, k);
    inc(p, k)
    ii = ii div wl
    inc(i)
  if (p and 1) != 0:
    var i: cint = 0
    while i < nd:
      var wl: cint = l.outerGeom[i]
      if i == l.innerCbDir: coords[i] = (coords[i] + l.innerCb * ib) mod wl
      inc(coords[i])
      if coords[i] >= wl:
        coords[i] = 0
        if i == l.innerCbDir: coords[i] = (coords[i] + l.innerCb * ib) mod wl
      else:
        if i == l.innerCbDir: coords[i] = (coords[i] + l.innerCb * ib) mod wl
        break
      inc(i)
  i=0
  while i < nd:
    inc(coords[i], cr[i])
    inc(i)
  var li2: LayoutIndexQ
  layoutIndexQ(l, addr(li2), coords)
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
