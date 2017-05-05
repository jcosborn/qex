import base
#import globals
#import basicOps
#import matrixConcept
import macros
#import metaUtils

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

template assignIadd(x,y:untyped):untyped = iadd(x,y)
template negIadd(x,y:untyped):untyped = isub(x,y)
template iaddIadd(x,y:untyped):untyped = iadd(x,y)
template isubIadd(x,y:untyped):untyped = isub(x,y)
template makeMap1(op:untyped):untyped {.dirty.} =
  template `op SV`*(rr:untyped; xx:untyped):untyped =
    subst(r,rr,x,xx,i,_):
      op(r, x[0])
      forO i, 1, <x.len:
        `op Iadd`(r, x[i])
  template `op SM`*(rr:typed; xx:typed):untyped =
    subst(r,rr,x,xx,i,_):
      assert(x.nrows == x.ncols)
      op(r, x[0,0])
      forO i, 1, <x.nrows:
        `op Iadd`(r, x[i,i])
  template `op VS`*(rr:typed; xx:typed):untyped =
    subst(r,rr,x,xx,tx,_,i,_):
      load(tx, x)
      forO i, 0, <r.len:
        op(r[i], tx)
  template `op VV`*(rr:untyped; xx:untyped):untyped =
    mixin op
    subst(r,rr,x,xx,i,_):
      assert(r.len == x.len)
      forO i, 0, <r.len:
        op(r[i], x[i])
  template `op MS`*(rr:typed; xx:typed):untyped =
    subst(r,rr,x,xx,tx,_,i,_,j,_):
      #assert(r.nrows == r.ncols)
      load(tx, x)
      cfor i, 0, <r.nrows:
        cfor j, 0, <r.ncols:
          if i == j:
            op(r[i,j], tx)
          else:
            op(r[i,j], 0)
  template `op MV`*(rr:typed; xx:typed):untyped =
    subst(r,rr,x,xx,i,_,j,_):
      assert(r.nrows == x.len)
      assert(r.ncols == x.len)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], x[i])
          else:
            op(r[i,j], 0)
  template `op MM`*(rr:untyped; xx:untyped):untyped =
    mixin op
    subst(r,rr,x,xx,i,_,j,_):
      assert(r.nrows == x.nrows)
      assert(r.ncols == x.ncols)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          op(r[i,j], x[i,j])

makeMap1(assign)
makeMap1(neg)
makeMap1(iadd)
makeMap1(isub)

template makeMap2(op:untyped):untyped {.dirty.} =
  template `op VVS`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,ty,_,i,_):
      assert(r.len == x.len)
      load(ty, y)
      forO i, 0, <r.len:
        op(r[i], x[i], ty)
  template `op VSV`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,tx,_,i,_):
      assert(r.len == y.len)
      load(tx, x)
      forO i, 0, <r.len:
        op(r[i], tx, y[i])
  template `op VVV`*(rr:typed; xx,yy:typed):untyped =
    mixin op
    subst(r,rr,x,xx,y,yy,i,_):
      assert(r.len == y.len)
      assert(r.len == x.len)
      forO i, 0, <r.len:
        op(r[i], x[i], y[i])
  template `op MSS`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,tx,_,ty,_,i,_,j,_):
      assert(r.nrows == r.ncols)
      load(tx, x)
      load(ty, y)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], tx, ty)
          else:
            op(r[i,j], 0, 0)
  #proc `op MVS`*(r:Vany; x,y:any) {.inline.} =
  template `op MVS`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,ty,_,i,_,j,_):
      assert(r.nrows == x.len)
      assert(r.ncols == x.len)
      load(ty, y)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], x[i], ty)
          else:
            op(r[i,j], 0, 0)
  template `op MSV`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,tx,_,i,_,j,_):
      assert(r.nrows == y.len)
      assert(r.ncols == y.len)
      load(tx, x)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], tx, y[i])
          else:
            op(r[i,j], 0, 0)
  template `op MVV`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,i,_,j,_):
      assert(r.nrows == x.len)
      assert(r.ncols == x.len)
      assert(x.len == y.len)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], x[i], y[i])
          else:
            op(r[i,j], 0, 0)
  template `op MMS`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,ty,_,i,_,j,_):
      mixin op
      assert(r.nrows == r.ncols)
      assert(r.nrows == x.nrows)
      assert(r.ncols == x.ncols)
      load(ty, y)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], x[i,j], ty)
          else:
            op(r[i,j], x[i,j], 0)
  template `op MSM`*(rr:typed; xx,yy:typed):untyped =
    mixin op
    subst(r,rr,x,xx,y,yy,tx,_,i,_,j,_):
      assert(r.nrows == r.ncols)
      assert(r.nrows == y.nrows)
      assert(r.ncols == y.ncols)
      load(tx, x)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], tx, y[i,j])
          else:
            op(r[i,j], 0, y[i,j])
  template `op MMV`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,i,_,j,_):
      assert(r.nrows == r.ncols)
      assert(r.nrows == x.nrows)
      assert(r.ncols == x.ncols)
      assert(r.nrows == y.len)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], x[i,j], y[i])
          else:
            op(r[i,j], x[i,j], 0)
  template `op MVM`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,i,_,j,_):
      assert(r.nrows == r.ncols)
      assert(r.nrows == x.len)
      assert(r.nrows == y.nrows)
      assert(r.ncols == y.ncols)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          if i == j:
            op(r[i,j], x[i], y[i,j])
          else:
            op(r[i,j], 0, y[i,j])
  template `op MMM`*(rr:typed; xx,yy:typed):untyped =
    mixin op
    subst(r,rr,x,xx,y,yy,i,_,j,_):
      assert(r.nrows == x.nrows)
      assert(r.ncols == x.ncols)
      assert(r.nrows == y.nrows)
      assert(r.ncols == y.ncols)
      forO i, 0, <r.nrows:
        forO j, 0, <r.ncols:
          op(r[i,j], x[i,j], y[i,j])

makeMap2(add)
makeMap2(sub)

template imulVS*(rr:untyped; xx:untyped):untyped =
  mixin imul
  subst(r,rr,x,xx,i,_,tx,_):
    load(tx, x)
    forO i, 0, <r.len:
      imul(r[i], tx)

template imulMS*(rr:untyped; xx:untyped):untyped =
  mixin imul
  subst(r,rr,x,xx,i,_,j,_,tx,_):
    load(tx, x)
    forO i, 0, <r.nrows:
      forO j, 0, <r.ncols:
        imul(r[i,j], tx)

template mulSVV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,i,_):
    mixin mul, imadd, assign
    assert(x.len == y.len)
    tmpvar(tr, r)
    mul(tr, x[0], y[0])
    forO i, 1, <x.len:
      imadd(tr, x[i], y[i])
    assign(r, tr)

template mulVVS*(r:typed; x,y:typed):untyped =
  mixin mul
  assert(r.len == x.len)
  load(ty, y)
  forO i, 0, <r.len:
    mul(r[i], x[i], ty)

template mulVSV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tx,_,i,_):
    mixin load, mul
    assert(r.len == y.len)
    load(tx, x)
    forO i, 0, <r.len:
      mul(r[i], tx, y[i])

template mulMMS*(rr:untyped; xx,yy:untyped):untyped =
  subst(r,rr,x,xx,y,yy,tx,_,i,_,j,_):
    #mixin mul
    assert(r.nrows == x.nrows)
    assert(r.ncols == x.ncols)
    load(ty, y)
    forO i, 0, <r.nrows:
      forO j, 0, <r.ncols:
        mul(r[i,j], x[i,j], ty)

template mulMSM*(rr:typed; xx,yy:typed):untyped =
  mixin mul
  subst(r,rr,x,xx,y,yy,tx,_,i,_,j,_):
    assert(r.nrows == y.nrows)
    assert(r.ncols == y.ncols)
    load(tx, x)
    forO i, 0, <r.nrows:
      forO j, 0, <r.ncols:
        mul(r[i,j], tx, y[i,j])

template mulVMV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,ty,_,i,_,j,_,ty0r,_,ty0i,_,tyjr,_,tyji,_):
    mixin nrows, ncols, mul, imadd, assign, load1
    assert(x.nrows == r.len)
    assert(x.ncols == y.len)
    when true:
    #when false:
      tmpvar(tr, r)
      #var tr{.noInit.}:array[r.len,type(load1(r[0]))]
      load(ty0r, y[0].re)
      forO i, 0, <x.nrows:
        mulCCR(tr[i], x[i,0], ty0r)
      load(ty0i, y[0].im)
      forO i, 0, <x.nrows:
        imaddCCI(tr[i], x[i,0], ty0i)
      forO j, 1, <x.ncols:
        load(tyjr, y[j].re)
        forO i, 0, <x.nrows:
          imaddCCR(tr[i], x[i,j], tyjr)
        load(tyji, y[j].im)
        forO i, 0, <x.nrows:
          imaddCCI(tr[i], x[i,j], tyji)
      assign(r, tr)
      #forO i, 0, <r.len: assign(r[i], tr[i])
    else:
      tmpvar(tr, r)
      block:
        load(ty, y[0])
        forO i, 0, <x.nrows:
          mul(tr[i], x[i,0], ty)
      forO j, 1, <x.ncols:
        load(ty, y[j])
        forO i, 0, <x.nrows:
          imadd(tr[i], x[i,j], ty)
      assign(r, tr)

template mulMMM*(rr:typed; xx,yy:typed):untyped =
  #[
  subst(r,rr,x,xx,y,yy,tr,_,i,_,j,_,k,_,txi0r,_,txi0i,_,txikr,_,txiki,_):
    mixin nrows, ncols, mul, imadd, assign, load1
    assert(r.nrows == x.nrows)
    assert(r.ncols == y.ncols)
    assert(x.ncols == y.nrows)
    var tr{.noInit.}:VectorArray[r.ncols,type(x[0,0]*y[0,0])]
    forO i, 0, <r.nrows:
      load(txi0r, x[i,0].re)
      forO j, 0, <r.ncols:
        mulCRC(tr[j], txi0r, y[0,j])
      load(txi0i, x[i,0].im)
      forO j, 0, <r.ncols:
        imaddCIC(tr[j], txi0i, y[0,j])
      forO k, 1, <x.ncols:
        load(txikr, x[i,k].re)
        forO j, 0, <r.ncols:
          imaddCRC(tr[j], txikr, y[k,j])
        load(txiki, x[i,k].im)
        forO j, 0, <r.ncols:
          imaddCIC(tr[j], txiki, y[k,j])
      forO j, 0, <r.ncols:
        assign(r[i,j], tr[j])
  ]#
  subst(r,rr,x,xx,y,yy,tr,_,i,_,j,_,k,_,txi0r,_,txi0i,_,txikr,_,txiki,_):
    assert(x.nrows == r.nrows)
    assert(x.ncols == y.nrows)
    assert(r.ncols == y.ncols)
    mixin mul, imadd
    forO i, 0, <r.nrows:
      var tr{.noInit.}:VectorArray[r.ncols,type(x[0,0]*y[0,0])]
      load(txi0, x[i,0])
      forO j, 0, <r.ncols:
        mul(tr[j], txi0, y[0,j])
      forO k, 1, <x.ncols:
        load(txik, x[i,k])
        forO j, 0, <r.ncols:
          imadd(tr[j], txik, y[k,j])
      forO j, 0, <r.ncols:
        assign(r[i,j], tr[j])


template imaddSVV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,i,_):
    mixin imadd, assign
    assert(x.len == y.len)
    load(tr, r)
    forO i, 0, <x.len:
      imadd(tr, x[i], y[i])
    assign(r, tr)

template imaddVMV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,ty,_,i,_,j,_,tyjr,_,tyji,_):
    mixin nrows, ncols, mul, imadd, assign, load1
    assert(x.nrows == r.len)
    assert(x.ncols == y.len)
    when true:
    #when false:
      load(tr, r)
      #var tr{.noInit.}:array[r.len,type(load1(r[0]))]
      #forO i, 0, <r.len: assign(tr[i], r[i])
      forO j, 0, <x.ncols:
        load(tyjr, y[j].re)
        forO i, 0, <x.nrows:
          imaddCCR(tr[i], x[i,j], tyjr)
        load(tyji, y[j].im)
        forO i, 0, <x.nrows:
          imaddCCI(tr[i], x[i,j], tyji)
      assign(r, tr)
      #forO i, 0, <r.len: assign(r[i], tr[i])
    else:
      load(tr, r)
      forO j, 0, <x.ncols:
        load(tyr, asReal(y[j].re))
        forO i, 0, <x.nrows:
          imadd(tr[i], x[i,j], tyr)
        load(tyi, asImag(y[j].im))
        forO i, 0, <x.nrows:
          imadd(tr[i], x[i,j], tyi)
      assign(r, tr)

template imaddMMM*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,ty,_,i,_,j,_,k,_,tyjr,_,tyji,_):
    mixin nrows, ncols, mul, imadd, assign, load1
    assert(r.nrows == x.nrows)
    assert(r.ncols == y.ncols)
    assert(x.ncols == y.nrows)
    #when true:
    when false:
      load(tr, r)
      forO i, 0, <r.nrows:
        forO k, 0, <x.ncols:
          load(txikr, x[i,k].re)
          forO j, 0, <r.ncols:
            imaddCRC(tr[i,j], txikr, y[k,j])
          load(txiki, x[i,k].im)
          forO j, 0, <r.ncols:
            imaddCIC(tr[i,j], txiki, y[k,j])
      assign(r, tr)
    else:
      forO i, 0, <r.nrows:
        var tr{.noInit.}:VectorArray[r.ncols,type(x[0,0]*y[0,0])]
        forO j, 0, <r.ncols:
          assign(tr[j], r[i,j])
        forO k, 0, <x.ncols:
          load(txikr, x[i,k].re)
          forO j, 0, <r.ncols:
            imaddCRC(tr[j], txikr, y[k,j])
          load(txiki, x[i,k].im)
          forO j, 0, <r.ncols:
            imaddCIC(tr[j], txiki, y[k,j])
        forO j, 0, <r.ncols:
          assign(r[i,j], tr[j])

template imsubVSV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tx,_,i,_):
    mixin imsub
    assert(r.len == y.len)
    load(tx, x)
    forO i, 0, <r.len:
      imsub(r[i], x, y[i])

template imsubVMV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,ty,_,i,_,j,_):
    mixin imsub
    assert(x.nrows == r.len)
    assert(x.ncols == y.len)
    when true:
      load(tr, r)
      forO j, 0, <x.ncols:
        load(ty, y[j])
        forO i, 0, <x.nrows:
          imsub(tr[i], x[i,j], ty)
      assign(r, tr)
    else:
      load(tr, r)
      forO j, 0, <x.ncols:
        load(tyr, asReal(y[j].re))
        forO i, 0, <x.nrows:
          imsub(tr[i], x[i,j], tyr)
        load(tyi, asImag(y[j].im))
        forO i, 0, <x.nrows:
          imsub(tr[i], x[i,j], tyi)
      assign(r, tr)

template imsubMMM*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,ty,_,i,_,j,_,k,_,txikr,_,txiki,_):
    mixin nrows, ncols, imsubCRC, imsubCIC, assign, load1
    assert(r.nrows == x.nrows)
    assert(r.ncols == y.ncols)
    assert(x.ncols == y.nrows)
    load(tr, r)
    forO i, 0, <r.nrows:
      forO k, 0, <x.ncols:
        load(txikr, x[i,k].re)
        forO j, 0, <r.ncols:
          imsubCRC(tr[i,j], txikr, y[k,j])
        load(txiki, x[i,k].im)
        forO j, 0, <r.ncols:
          imsubCIC(tr[i,j], txiki, y[k,j])
    assign(r, tr)

template msubVSVV*(rr:typed; xx,yy,zz:typed):untyped =
  subst(r,rr,x,xx,y,yy,z,zz,i,_):
    mixin msub
    assert(r.len == y.len)
    assert(r.len == z.len)
    forO i, 0, <r.len:
      msub(r[i], x, y[i], z[i])
