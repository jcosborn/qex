import basicOps
#import matrixConcept
import macros
import metaUtils

template cfor*(i,r0,r1,b:untyped):untyped =
  block:
    var i = r0
    while i <= r1:
      b
      inc(i)
#template forO*(i,r0,r1,b:untyped):untyped = #cfor(i,r0,r1,b)
#  var i:int
#  for ii{.gensym.} in r0..r1:
#    i = ii
#    b
macro forO*(i,r0,r1,b:untyped):auto =
  #echo b.repr
  result = quote do:
    for `i` in `r0`..`r1`:
      `b`
#template forO*(i,r0,r1,b:untyped):untyped = forStatic(i,r0,r1,b)

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

proc mulMMS*(r:any; x,y:any) {.inline.} =
  #mixin mul
  assert(r.nrows == x.nrows)
  assert(r.ncols == x.ncols)
  #load(ty, y)
  forO i, 0, <r.nrows:
    forO j, 0, <r.ncols:
      #echo isComplex(r[i,j])
      mul(r[i,j], x[i,j], y)

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
        block:
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
