import macros
import base/globals
import base/basicOps
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
  template `op VS`*(r: typed; xx: typed) =
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
      #static: echo "opVS: ", x.type, " ", xx.type, " ", astToStr(xx)
      #forO i, 0, r.len.pred:
      for i in fOpt(0,r.len.pred):
        op(r[i], x)
  template `op VV`*(r: typed; xx: typed) =
    #static: echo instantiationInfo(0)
    #static: echo instantiationInfo(1)
    #static: echo instantiationInfo(2)
    #static: echo instantiationInfo(3)
    #static: echo instantiationInfo(-1)
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
      #static: echo "opVV: ", astToStr(op), " ", x.type, " ", xx.type, " ", astToStr(xx)
      assert(r.len == x.len)
      #forO i, 0, r.len.pred:
      for i in fOpt(0,r.len.pred):
        op(r[i], x[i])
        #let xi = x[i]   ## FIXME
        #op(r[i], xi)
  template `op MS`*(r: typed; xx: typed) =
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
      #forO i, 0, r.nrows.pred:
      for i in fOpt(0,r.nrows.pred):
        #forO j, 0, r.ncols.pred:
        for j in fOpt(0,r.ncols.pred):
          if i == j:
            op(r[i,j], x)
          else:
            op(r[i,j], 0)
  template `op MV`*(r:typed; xx:typed) =
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
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
  template `op MM`*(r:typed; xx:typed) =
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
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
  template `op VVS`*(r:typed; xx,yy:typed) =
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
      let yp = getPtr yy; template y:untyped = yp[]
      assert(r.len == x.len)
      #forO i, 0, r.len.pred:
      for i in fOpt(0,r.len.pred):
        op(r[i], x[i], y)
  template `op VSV`*(r:typed; xx,yy:typed) =
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
      let yp = getPtr yy; template y:untyped = yp[]
      assert(r.len == y.len)
      #forO i, 0, r.len.pred:
      for i in fOpt(0,r.len.pred):
        op(r[i], x, y[i])
  template `op VVV`*(r: typed; xx,yy: typed) =
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
      let yp = getPtr yy; template y:untyped = yp[]
      assert(r.len == y.len)
      assert(r.len == x.len)
      #forO i, 0, r.len.pred:
      for i in fOpt(0,r.len.pred):
        op(r[i], x[i], y[i])
  template `op MSS`*(r:typed; xx,yy:typed) =
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
      let yp = getPtr yy; template y:untyped = yp[]
      assert(r.nrows == r.ncols)
      #forO i, 0, r.nrows.pred:
      for i in fOpt(0,r.nrows.pred):
        #forO j, 0, r.ncols.pred:
        for j in fOpt(0,r.ncols.pred):
          if i == j:
            op(r[i,j], x, y)
          else:
            op(r[i,j], 0, 0)
  template `op MVS`*(r:typed; xx,yy:typed) =
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
      let yp = getPtr yy; template y:untyped = yp[]
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
  template `op MSV`*(r:typed; xx,yy:typed) =
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
      let yp = getPtr yy; template y:untyped = yp[]
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
  template `op MVV`*(r:typed; xx,yy:typed) =
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
      let yp = getPtr yy; template y:untyped = yp[]
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
  template `op MMS`*(r:typed; xx,yy:typed) =
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
      let yp = getPtr yy; template y:untyped = yp[]
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
  template `op MSM`*(r:typed; xx,yy:typed) =
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
      let yp = getPtr yy; template y:untyped = yp[]
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
  template `op MMV`*(r:typed; xx,yy:typed) =
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
      let yp = getPtr yy; template y:untyped = yp[]
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
  template `op MVM`*(r:typed; xx,yy:typed) =
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
      let yp = getPtr yy; template y:untyped = yp[]
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
  template `op MMM`*(r:typed; xx,yy:typed) =
    mixin op
    block:
      let xp = getPtr xx; template x:untyped = xp[]
      let yp = getPtr yy; template y:untyped = yp[]
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
    let xp = getPtr xx; template x:untyped = xp[]
    #forO i, 0, r.len.pred:
    for i in fOpt(0,r.len.pred):
      imul(r[i], x)

template imulMS*(r: typed; xx: typed) =
  mixin imul
  block:
    let xp = getPtr xx; template x:untyped = xp[]
    #forO i, 0, r.nrows.pred:
    for i in fOpt(0,r.nrows.pred):
      #forO j, 0, r.ncols.pred:
      for j in fOpt(0,r.ncols.pred):
        imul(r[i,j], x)

template mulSVV*(r:typed; xx,yy:typed) =
  mixin mul, imadd #, assign
  block:
    let xp = getPtr xx; template x:untyped = xp[]
    let yp = getPtr yy; template y:untyped = yp[]
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
    let xp = getPtr xx; template x:untyped = xp[]
    let yp = getPtr yy; template y:untyped = yp[]
    assert(r.len == x.len)
    #forO i, 0, r.len.pred:
    for i in fOpt(0,r.len.pred):
      mul(r[i], x[i], y)

template mulVSV*(r:typed; xx,yy:typed) =
  mixin mul
  block:
    let xp = getPtr xx; template x:untyped = xp[]
    let yp = getPtr yy; template y:untyped = yp[]
    assert(r.len == y[].len)
    #forO i, 0, r.len.pred:
    for i in fOpt(0,r.len.pred):
      mul(r[i], x, y[i])

template mulMMS*(r:typed; xx,yy:typed) =
  mixin mul
  block:
    let xp = getPtr xx; template x:untyped = xp[]
    let yp = getPtr yy; template y:untyped = yp[]
    assert(r.nrows == x.nrows)
    assert(r.ncols == x.ncols)
    #forO i, 0, r.nrows.pred:
    for i in fOpt(0,r.nrows.pred):
      #forO j, 0, r.ncols.pred:
      for j in fOpt(0,r.ncols.pred):
        mul(r[i,j], x[i,j], y)

template mulMSM*(r:typed; xx,yy:typed) =
  mixin mul
  block:
    let xp = getPtr xx; template x:untyped = xp[]
    let yp = getPtr yy; template y:untyped = yp[]
    assert(r.nrows == y.nrows)
    assert(r.ncols == y.ncols)
    #forO i, 0, r.nrows.pred:
    for i in fOpt(0,r.nrows.pred):
      #forO j, 0, r.ncols.pred:
      for j in fOpt(0,r.ncols.pred):
        mul(r[i,j], x, y[i,j])

template mulVMV*(r: typed; xx,yy: typed) =
  mixin nrows, ncols, mul, imadd
  block:
    let xp = getPtr xx; template x:untyped = xp[]
    let yp = getPtr yy; template y:untyped = yp[]
    assert(x.nrows == r.len)
    assert(x.ncols == y.len)
    #forO i, 0, x.nrows.pred:
    for i in fOpt(0,x.nrows.pred):
      mul(r[i], x[i,0], y[0])
      #forO j, 1, x.ncols.pred:
      for j in fOpt(1,x.ncols.pred):
        imadd(r[i], x[i,j], y[j])

template mulMMM*(r: typed; xx,yy: typed) =
  mixin mul, imadd
  block:
    let xp = getPtr xx; template x:untyped = xp[]
    let yp = getPtr yy; template y:untyped = yp[]
    assert(x.nrows == r.nrows)
    assert(x.ncols == y.nrows)
    assert(r.ncols == y.ncols)
    forO i, 0, r.nrows.pred:
      forO j, 0, r.ncols.pred:
        mul(r[i,j], x[i,0], y[0,j])
        forO k, 1, y.nrows.pred:
          imadd(r[i,j], x[i,k], y[k,j])

template imaddSVV*(r:typed; xx,yy:typed) =
  mixin imadd, assign
  block:
    let xp = getPtr xx; template x:untyped = xp[]
    let yp = getPtr yy; template y:untyped = yp[]
    assert(x.len == y.len)
    forO i, 0, x.len.pred:
      imadd(r, x[i], y[i])

template imaddVSV*(r: typed; xx,yy: typed) =
  mixin imadd
  block:
    let xp = getPtr xx; template x:untyped = xp[]
    let yp = getPtr yy; template y:untyped = yp[]
    assert(r.len == y.len)
    forO i, 0, r.len.pred:
      imadd(r[i], x, y[i])

template imaddVMV*(r: typed; xx,yy: typed) =
  mixin nrows, ncols, imadd
  block:
    let xp = getPtr xx; template x:untyped = xp[]
    let yp = getPtr yy; template y:untyped = yp[]
    assert(x.nrows == r.len)
    assert(x.ncols == y.len)
    for j in fOpt(0,x.ncols.pred):
      for i in fOpt(0,x.nrows.pred):
        imadd(r[i], x[i,j], y[j])

template imaddMMM*(r:typed; xx,yy:typed) =
  mixin nrows, ncols, imadd
  block:
    let xp = getPtr xx; template x:untyped = xp[]
    let yp = getPtr yy; template y:untyped = yp[]
    assert(r.nrows == x.nrows)
    assert(r.ncols == y.ncols)
    assert(x.ncols == y.nrows)
    for i in fOpt(0,r.nrows.pred):
      for j in fOpt(0,r.ncols.pred):
        for k in fOpt(0,x.ncols.pred):
          imadd(r[i,j], x[i,k], y[k,j])

template imsubVSV*(r:typed; xx,yy:typed) =
  mixin imsub
  block:
    let xp = getPtr xx; template x:untyped = xp[]
    let yp = getPtr yy; template y:untyped = yp[]
    assert(r.len == y.len)
    for i in fOpt(0,r.len.pred):
      imsub(r[i], x, y[i])

template imsubVMV*(r:typed; xx,yy:typed) =
  mixin imsub
  block:
    let xp = getPtr xx; template x:untyped = xp[]
    let yp = getPtr yy; template y:untyped = yp[]
    assert(x.nrows == r.len)
    assert(x.ncols == y.len)
    for j in fOpt(0,x.ncols.pred):
      for i in fOpt(0,x.nrows.pred):
        imsub(r[i], x[i,j], y[j])

template imsubMMM*(r:typed; xx,yy:typed) =
  mixin nrows, ncols, imsub
  block:
    let xp = getPtr xx; template x:untyped = xp[]
    let yp = getPtr yy; template y:untyped = yp[]
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
    let xp = getPtr xx; template x:untyped = xp[]
    let yp = getPtr yy; template y:untyped = yp[]
    let zp = getPtr zz; template z:untyped = zp[]
    assert(r.len == y.len)
    assert(r.len == z.len)
    for i in fOpt(0,r.len.pred):
      msub(r[i], x, y[i], z[i])
