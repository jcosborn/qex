import macros
import base/globals
import base/basicOps
import base/wrapperTypes
import base/metaUtils
#import globals
#import basicOps
#import matrixConcept

template cfor*(i,r0,r1,b:untyped):untyped =
  block:
    var i = r0
    while i <= r1:
      b
      inc i
#template forO*(i,r0,r1,b:untyped):untyped = #cfor(i,r0,r1,b)
#  var i:int
#  for ii{.gensym.} in r0..r1:
#    i = ii
#    b
macro forN*(i,r0,r1,b:untyped):auto =
  #echo b.repr
  result = quote do:
    #for `i` in `r0`..`r1`:
    for `i` in countup(`r0`,`r1`):
      `b`
when staticUnroll:
  template forO*(i,r0,r1,b: untyped): untyped {.dirty.} =
    bind forStatic
    forStatic(i,r0,r1,b)
else:
  template forO*(i,r0,r1,b:untyped):untyped {.dirty.} =
    forN(i,r0,r1,b)
  #template forO*(i,r0,r1,b:untyped):untyped = cfor(i,r0,r1,b)

macro fOpt(stmt: ForLoopStmt): untyped =
  let expr = stmt[0]
  let iter = stmt[1]
  let a = iter[1]
  let b = iter[2]
  let body = stmt[2]
  result = quote do:
    forO `expr`, `a`, `b`:
      `body`

#template assignIadd(x,y:typed) = iadd(x,y)
#template negIadd(x,y:typed) = isub(x,y)
#template iaddIadd(x,y:typed) = iadd(x,y)
#template isubIadd(x,y:typed) = isub(x,y)
template makeMap1(op:untyped) =
  getOptimPragmas()
  template `op VS`*(rr: typed; xx: typed) =
    mixin op
    block:
      let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
      let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
      #let x {.byptr.} = xx
      #static: echo "opVS: ", x.type, " ", xx.type, " ", astToStr(xx)
      #forO i, 0, r.len.pred:
      for i in fOpt(0,r.len.pred):
        op(r[i], x)
  proc `op VV`*(r: var auto; x: auto) {.alwaysInline.} =
  #template `op VV`*(rr: typed; xx: typed) =
  #  let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
  #  let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
    mixin op
    assert(r.len == x.len)
    #forO i, 0, r.len.pred:
    for i in fOpt(0,r.len.pred):
      op(r[i], x[i])
  #[
  template `op MS`*(rr: typed; xx: typed) =
    mixin op
    block:
      let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
      let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
      #let x {.byptr.} = xx
      #forO i, 0, r.nrows.pred:
      for i in fOpt(0,r.nrows.pred):
        #forO j, 0, r.ncols.pred:
        for j in fOpt(0,r.ncols.pred):
          if i == j:
            op(r[i,j], x)
          else:
            op(r[i,j], 0)
  ]#
  proc `op MS`*(r: var auto; x: auto) {.alwaysInline.} =
    mixin op
    when astToStr(op) == "iadd" or astToStr(op) == "isub":
      for i in fOpt(0,min(r.nrows.pred,r.ncols.pred)):
        op(r[i,i], x)
    else:
      #forO i, 0, r.nrows.pred:
      for i in fOpt(0,r.nrows.pred):
        #forO j, 0, r.ncols.pred:
        for j in fOpt(0,r.ncols.pred):
          if i == j:
            op(r[i,j], x)
          else:
            op(r[i,j], 0)
  template `op MV`*(rr:typed; xx:typed) =
    mixin op
    block:
      let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
      let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
      #let x {.byptr.} = xx
      assert(r.nrows == x.len)
      assert(r.ncols == x.len)
      #forO i, 0, r.nrows.pred:
      for i in fOpt(0,r.nrows.pred):
        #forO j, 0, r.ncols.pred:
        for j in fOpt(0,r.ncols.pred):
          if i == j:
            op(r[i,j], x[i])
          else:
            op(r[i,j], 0)
  #[
  template `op MM`*(rr:typed; xx:typed) =
    mixin op
    block:
      let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
      let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
      #let x {.byptr.} = xx
      assert(r.nrows == x.nrows)
      assert(r.ncols == x.ncols)
      #forO i, 0, r.nrows.pred:
      for i in fOpt(0,r.nrows.pred):
        #forO j, 0, r.ncols.pred:
        for j in fOpt(0,r.ncols.pred):
          op(r[i,j], x[i,j])
  ]#
  proc `op MM`*(r: var auto; x: auto) {.alwaysInline.} =
    mixin op
    assert(r.nrows == x.nrows)
    assert(r.ncols == x.ncols)
    #forO i, 0, r.nrows.pred:
    for i in fOpt(0,r.nrows.pred):
      #forO j, 0, r.ncols.pred:
      for j in fOpt(0,r.ncols.pred):
        op(r[i,j], x[i,j])

makeMap1(assign)
makeMap1(neg)
makeMap1(iadd)
makeMap1(isub)

template makeMap2(op:untyped):untyped {.dirty.} =
  getOptimPragmas()
  template `op VVS`*(rr:typed; xx,yy:typed) =
    mixin op
    block:
      let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
      let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
      let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
      assert(r.len == x.len)
      #forO i, 0, r.len.pred:
      for i in fOpt(0,r.len.pred):
        op(r[i], x[i], y)
  template `op VSV`*(rr:typed; xx,yy:typed) =
    mixin op
    block:
      let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
      let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
      let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
      assert(r.len == y.len)
      #forO i, 0, r.len.pred:
      for i in fOpt(0,r.len.pred):
        op(r[i], x, y[i])
  proc `op VVV`*(r: var auto; x,y: auto) {.alwaysInline.} =
    mixin op
    assert(r.len == y.len)
    assert(r.len == x.len)
    #forO i, 0, r.len.pred:
    for i in fOpt(0,r.len.pred):
      op(r[i], x[i], y[i])
  template `op MSS`*(rr:typed; xx,yy:typed) =
    mixin op
    block:
      let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
      let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
      let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
      assert(r.nrows == r.ncols)
      #forO i, 0, r.nrows.pred:
      for i in fOpt(0,r.nrows.pred):
        #forO j, 0, r.ncols.pred:
        for j in fOpt(0,r.ncols.pred):
          if i == j:
            op(r[i,j], x, y)
          else:
            op(r[i,j], 0, 0)
  template `op MVS`*(rr:typed; xx,yy:typed) =
    mixin op
    block:
      let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
      let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
      let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
      assert(r.nrows == x.len)
      assert(r.ncols == x.len)
      #forO i, 0, r.nrows.pred:
      for i in fOpt(0,r.nrows.pred):
        #forO j, 0, r.ncols.pred:
        for j in fOpt(0,r.ncols.pred):
          if i == j:
            op(r[i,j], x[i], y)
          else:
            op(r[i,j], 0, 0)
  template `op MSV`*(rr:typed; xx,yy:typed) =
    mixin op
    block:
      let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
      let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
      let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
      assert(r.nrows == y.len)
      assert(r.ncols == y.len)
      #forO i, 0, r.nrows.pred:
      for i in fOpt(0,r.nrows.pred):
        #forO j, 0, r.ncols.pred:
        for j in fOpt(0,r.ncols.pred):
          if i == j:
            op(r[i,j], x, y[i])
          else:
            op(r[i,j], 0, 0)
  template `op MVV`*(rr:typed; xx,yy:typed) =
    mixin op
    block:
      let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
      let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
      let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
      assert(r.nrows == x.len)
      assert(r.ncols == x.len)
      assert(x.len == y.len)
      #forO i, 0, r.nrows.pred:
      for i in fOpt(0,r.nrows.pred):
        #forO j, 0, r.ncols.pred:
        for j in fOpt(0,r.ncols.pred):
          if i == j:
            op(r[i,j], x[i], y[i])
          else:
            op(r[i,j], 0, 0)
  #[
  template `op MMS`*(rr:typed; xx,yy:typed) =
    mixin op
    block:
      let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
      let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
      let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
      assert(r.nrows == r.ncols)
      assert(r.nrows == x.nrows)
      assert(r.ncols == x.ncols)
      #forO i, 0, r.nrows.pred:
      for i in fOpt(0,r.nrows.pred):
        #forO j, 0, r.ncols.pred:
        for j in fOpt(0,r.ncols.pred):
          if i == j:
            op(r[i,j], x[i,j], y)
          else:
            op(r[i,j], x[i,j], 0)
  ]#
  proc `op MMS`*(r: var auto; x,y: auto) {.alwaysInline.} =
    mixin op
    assert(r.nrows == r.ncols)
    assert(r.nrows == x.nrows)
    assert(r.ncols == x.ncols)
    #forO i, 0, r.nrows.pred:
    for i in fOpt(0,r.nrows.pred):
      #forO j, 0, r.ncols.pred:
      for j in fOpt(0,r.ncols.pred):
        if i == j:
          op(r[i,j], x[i,j], y)
        else:
          op(r[i,j], x[i,j], 0)
  #[
  template `op MSM`*(rr:typed; xx,yy:typed) =
    mixin op
    block:
      let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
      let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
      let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
      assert(r.nrows == r.ncols)
      assert(r.nrows == y.nrows)
      assert(r.ncols == y.ncols)
      #forO i, 0, r.nrows.pred:
      for i in fOpt(0,r.nrows.pred):
        #forO j, 0, r.ncols.pred:
        for j in fOpt(0,r.ncols.pred):
          if i == j:
            op(r[i,j], x, y[i,j])
          else:
            op(r[i,j], 0, y[i,j])
  ]#
  proc `op MSM`*(r: var auto; x,y: auto) {.alwaysInline.} =
    mixin op
    assert(r.nrows == r.ncols)
    assert(r.nrows == y.nrows)
    assert(r.ncols == y.ncols)
    #forO i, 0, r.nrows.pred:
    for i in fOpt(0,r.nrows.pred):
      #forO j, 0, r.ncols.pred:
      for j in fOpt(0,r.ncols.pred):
        if i == j:
          op(r[i,j], x, y[i,j])
        else:
          op(r[i,j], 0, y[i,j])
  template `op MMV`*(rr:typed; xx,yy:typed) =
    mixin op
    block:
      let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
      let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
      let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
      assert(r.nrows == r.ncols)
      assert(r.nrows == x.nrows)
      assert(r.ncols == x.ncols)
      assert(r.nrows == y.len)
      #forO i, 0, r.nrows.pred:
      for i in fOpt(0,r.nrows.pred):
        #forO j, 0, r.ncols.pred:
        for j in fOpt(0,r.ncols.pred):
          if i == j:
            op(r[i,j], x[i,j], y[i])
          else:
            op(r[i,j], x[i,j], 0)
  template `op MVM`*(rr:typed; xx,yy:typed) =
    mixin op
    block:
      let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
      let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
      let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
      assert(r.nrows == r.ncols)
      assert(r.nrows == x.len)
      assert(r.nrows == y.nrows)
      assert(r.ncols == y.ncols)
      #forO i, 0, r.nrows.pred:
      for i in fOpt(0,r.nrows.pred):
        #forO j, 0, r.ncols.pred:
        for j in fOpt(0,r.ncols.pred):
          if i == j:
            op(r[i,j], x[i], y[i,j])
          else:
            op(r[i,j], 0, y[i,j])
  #[
  template `op MMM`*(rr:typed; xx,yy:typed) =
    mixin op
    block:
      let rp = getPtr rr; template r:untyped {.gensym.} = rp[]
      let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
      let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
      assert(r.nrows == x.nrows)
      assert(r.ncols == x.ncols)
      assert(r.nrows == y.nrows)
      assert(r.ncols == y.ncols)
      #forO i, 0, r.nrows.pred:
      for i in fOpt(0,r.nrows.pred):
        #forO j, 0, r.ncols.pred:
        for j in fOpt(0,r.ncols.pred):
          op(r[i,j], x[i,j], y[i,j])
  ]#
  proc `op MMM`*(r: var auto; x,y: auto) {.alwaysInline.} =
    mixin op
    assert(r.nrows == x.nrows)
    assert(r.ncols == x.ncols)
    assert(r.nrows == y.nrows)
    assert(r.ncols == y.ncols)
    #forO i, 0, r.nrows.pred:
    for i in fOpt(0,r.nrows.pred):
      #forO j, 0, r.ncols.pred:
      for j in fOpt(0,r.ncols.pred):
        op(r[i,j], x[i,j], y[i,j])

makeMap2(add)
makeMap2(sub)

template imulVS*(r:typed; xx:typed) =
  mixin imul
  block:
    let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
    #forO i, 0, r.len.pred:
    for i in fOpt(0,r.len.pred):
      imul(r[i], x)

#[
template imulMS*(r: typed; xx: typed) =
  mixin imul
  block:
    let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
    #forO i, 0, r.nrows.pred:
    for i in fOpt(0,r.nrows.pred):
      #forO j, 0, r.ncols.pred:
      for j in fOpt(0,r.ncols.pred):
        imul(r[i,j], x)
]#
proc imulMS*(r: var auto; x: auto) {.alwaysInline.}=
  mixin imul
  #forO i, 0, r.nrows.pred:
  for i in fOpt(0,r.nrows.pred):
    #forO j, 0, r.ncols.pred:
    for j in fOpt(0,r.ncols.pred):
      imul(r[i,j], x)

template mulSVV*(r:typed; xx,yy:typed) =
  mixin mul, imadd #, assign
  block:
    let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
    let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
    assert(x.len == y.len)
    #tmpvar(tr, r)
    mul(r, x[0], y[0])
    #forO i, 1, x.len.pred:
    for i in fOpt(1,x.len.pred):
      imadd(r, x[i], y[i])
    #assign(r, tr)

template mulVVS*(r:typed; xx,yy:typed) =
  mixin mul
  block:
    let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
    let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
    assert(r.len == x.len)
    #forO i, 0, r.len.pred:
    for i in fOpt(0,r.len.pred):
      mul(r[i], x[i], y)

proc mulVSV*(r: var auto; x,y: auto) {.alwaysInline.} =
  mixin mul
  assert(r.len == y[].len)
  #forO i, 0, r.len.pred:
  for i in fOpt(0,r.len.pred):
    mul(r[i], x, y[i])

#[
template mulVSV*(rr:typed; xx,yy:typed):untyped =
  #subst(r,rr,x,xx,y,yy,tx,_,i,_):
  subst(r,rr,x,xx,y,yy,i,_):
    mixin load1, mul
    assert(r.len == y.len)
    #load(tx, x)
    let txz = load1(x)
    forO i, 0, r.len.pred:
      mul(r[i], txz, y[i])
template mulVSVU*(r: typed; x,y: typed): untyped =
  mixin mul, `:=`
  assert(r.len == y.len)
  #forO i, 0, r.len.pred:
  for i in fOpt(0,r.len.pred):
    #mul(r[i], x, y[i])
    r[i] := mul(x, y[i])
template mulVSV*(r: typed; x,y: typed): untyped =
  # prepMultipleAccess(x)
  # prepLinearAccess(y)
  flattenCallArgs(mulVSVU, r, x, y)
]#

template mulMMS*(r:typed; xx,yy:typed) =
  mixin mul
  block:
    let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
    let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
    assert(r.nrows == x.nrows)
    assert(r.ncols == x.ncols)
    #forO i, 0, r.nrows.pred:
    for i in fOpt(0,r.nrows.pred):
      #forO j, 0, r.ncols.pred:
      for j in fOpt(0,r.ncols.pred):
        mul(r[i,j], x[i,j], y)

#[
template mulMSM*(r:typed; xx,yy:typed) =
  mixin mul
  block:
    let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
    let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
    assert(r.nrows == y.nrows)
    assert(r.ncols == y.ncols)
    #forO i, 0, r.nrows.pred:
    for i in fOpt(0,r.nrows.pred):
      #forO j, 0, r.ncols.pred:
      for j in fOpt(0,r.ncols.pred):
        mul(r[i,j], x, y[i,j])
]#
proc mulMSM*(r: var auto; x,y: auto) {.alwaysInline.} =
  mixin mul
  assert(r.nrows == y.nrows)
  assert(r.ncols == y.ncols)
  #forO i, 0, r.nrows.pred:
  for i in fOpt(0,r.nrows.pred):
    #forO j, 0, r.ncols.pred:
    for j in fOpt(0,r.ncols.pred):
      mul(r[i,j], x, y[i,j])

proc mulVMV*(r: var auto; x,y: auto) {.alwaysInline.} =
  mixin nrows, ncols, mul, imadd
  assert(x.nrows == r.len)
  assert(x.ncols == y.len)
  when true:
    forO i, 0, x.nrows.pred:
    #for i in fOpt(0,x.nrows.pred):
      #mul(r[i], x[i,0], y[0])
      var t {.noInit.}: evalType(r[i])
      mul(t, x[i,0], y[0])
      forO j, 1, x.ncols.pred:
      #for j in fOpt(1,x.ncols.pred):
        #imadd(r[i], x[i,j], y[j])
        imadd(t, x[i,j], y[j])
      r[i] := t
  else:
    for i in fOpt(0,x.nrows.pred):
      mul(r[i], x[i,0], y[0])
    for j in fOpt(1,x.ncols.pred):
      for i in fOpt(0,x.nrows.pred):
        imadd(r[i], x[i,j], y[j])

#[
template mulMMM*(r: typed; xx,yy: typed) =
  mixin mul, imadd
  block:
    let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
    let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
    assert(x.nrows == r.nrows)
    assert(x.ncols == y.nrows)
    assert(r.ncols == y.ncols)
    var t {.noInit.}: evalType(r)
    forO i, 0, r.nrows.pred:
      forO j, 0, r.ncols.pred:
        #mul(r[i,j], x[i,0], y[0,j])
        mul(t[i,j], x[i,0], y[0,j])
        forO k, 1, y.nrows.pred:
          #imadd(r[i,j], x[i,k], y[k,j])
          imadd(t[i,j], x[i,k], y[k,j])
    r := t
]#
proc mulMMM*(r: var auto; x,y: auto) {.alwaysInline.} =
  mixin mul, imadd
  assert(x.nrows == r.nrows)
  assert(x.ncols == y.nrows)
  assert(r.ncols == y.ncols)
  forO i, 0, r.nrows.pred:
    forO j, 0, r.ncols.pred:
      var t {.noInit.}: evalType(r[i,j])
      mul(t, x[i,0], y[0,j])
      forO k, 1, y.nrows.pred:
        imadd(t, x[i,k], y[k,j])
      r[i,j] = t

template imaddSVV*(r:typed; xx,yy:typed) =
  mixin imadd, assign
  block:
    let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
    let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
    assert(x.len == y.len)
    forO i, 0, x.len.pred:
      imadd(r, x[i], y[i])

template imaddVSV*(r: typed; xx,yy: typed) =
  mixin imadd
  block:
    let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
    let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
    assert(r.len == y.len)
    forO i, 0, r.len.pred:
      imadd(r[i], x, y[i])

#template imaddVMV*(r: typed; xx,yy: typed) =
proc imaddVMV*(r: var auto; x,y: auto) {.alwaysInline.} =
  mixin nrows, ncols, imadd
  assert(x.nrows == r.len)
  assert(x.ncols == y.len)
  for i in fOpt(0,x.nrows.pred):
    var t {.noinit.}: evalType(r[i])
    t := r[i]
    for j in fOpt(0,x.ncols.pred):
      imadd(t, x[i,j], y[j])
    r[i] := t

#[
template imaddMMM*(r:typed; xx,yy:typed) =
  mixin nrows, ncols, imadd
  block:
    let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
    let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
    assert(r.nrows == x.nrows)
    assert(r.ncols == y.ncols)
    assert(x.ncols == y.nrows)
    for i in fOpt(0,r.nrows.pred):
      for j in fOpt(0,r.ncols.pred):
        for k in fOpt(0,x.ncols.pred):
          imadd(r[i,j], x[i,k], y[k,j])
]#
proc imaddMMM*(r: var auto; x,y: auto) {.alwaysInline.} =
  mixin nrows, ncols, imadd
  assert(r.nrows == x.nrows)
  assert(r.ncols == y.ncols)
  assert(x.ncols == y.nrows)
  for i in fOpt(0,r.nrows.pred):
    for j in fOpt(0,r.ncols.pred):
      var t {.noInit.}: evalType(r[i,j])
      t := r[i,j]
      for k in fOpt(0,x.ncols.pred):
        #imadd(r[i,j], x[i,k], y[k,j])
        imadd(t, x[i,k], y[k,j])
      r[i,j] = t

template imsubVSV*(r:typed; xx,yy:typed) =
  mixin imsub
  block:
    let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
    let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
    assert(r.len == y.len)
    for i in fOpt(0,r.len.pred):
      imsub(r[i], x, y[i])

template imsubVMV*(r:typed; xx,yy:typed) =
  mixin imsub
  block:
    let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
    let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
    assert(x.nrows == r.len)
    assert(x.ncols == y.len)
    for j in fOpt(0,x.ncols.pred):
      for i in fOpt(0,x.nrows.pred):
        imsub(r[i], x[i,j], y[j])

template imsubMMM*(r:typed; xx,yy:typed) =
  mixin nrows, ncols, imsub
  block:
    let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
    let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
    assert(r.nrows == x.nrows)
    assert(r.ncols == y.ncols)
    assert(x.ncols == y.nrows)
    for i in fOpt(0,r.nrows.pred):
      for j in fOpt(0,r.ncols.pred):
        for k in fOpt(0,x.ncols.pred):
          imsub(r[i,j], x[i,k], y[k,j])

template msubVSVV*(r:typed; xx,yy,zz:typed) =
  mixin msub
  block:
    let xp = getPtr xx; template x:untyped {.gensym.} = xp[]
    let yp = getPtr yy; template y:untyped {.gensym.} = yp[]
    let zp = getPtr zz; template z:untyped {.gensym.} = zp[]
    assert(r.len == y.len)
    assert(r.len == z.len)
    for i in fOpt(0,r.len.pred):
      msub(r[i], x, y[i], z[i])
