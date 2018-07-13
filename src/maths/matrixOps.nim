import base
#import globals
#import basicOps
#import matrixConcept
import macros

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
      forO i, 1, x.len.pred:
        `op Iadd`(r, x[i])
  template `op SM`*(rr:typed; xx:typed):untyped =
    subst(r,rr,x,xx,i,_):
      assert(x.nrows == x.ncols)
      op(r, x[0,0])
      forO i, 1, x.nrows.pred:
        `op Iadd`(r, x[i,i])
  template `op VS`*(r: typed; x: typed): untyped =
    #subst(r,rr,x,xx,tx,_,i,_):
    #load(tx, x)
    let tx_opVS = x
    forO i, 0, r.len.pred:
      op(r[i], tx_opVS)
  #[
  template `op VV`*(r: typed; xx: typed): untyped =
    mixin op
    #[
    optimizeAst:
      subst(r,rr,x,xx,i,_):
        assert(r.len == x.len)
        forO i, 0, r.len.pred:
          op(r[i], x[i])
    ]#
    let x = xx
    assert(r.len == x.len)
    forO i, 0, r.len.pred:
      op(r[i], x[i])
  ]#
  template `op VVU`*(r: typed; x: typed): untyped =
    mixin op
    assert(r.len == x.len)
    forO i, 0, r.len.pred:
      op(r[i], x[i])
  template `op VV`*(r: typed; x: typed): untyped =
    flattenCallArgs(`op VVU`, r, x)
  template `op MS`*(r: typed; x: typed): untyped =
    #subst(r,rr,x,xx,tx,_,i,_,j,_):
    #assert(r.nrows == r.ncols)
    #load(tx, x)
    let tx_opMS = x
    forO i, 0, r.nrows.pred:
      forO j, 0, r.ncols.pred:
        if i == j:
          op(r[i,j], tx_opMS)
        else:
          op(r[i,j], 0)
  template `op MV`*(rr:typed; xx:typed):untyped =
    subst(r,rr,x,xx,i,_,j,_):
      assert(r.nrows == x.len)
      assert(r.ncols == x.len)
      forO i, 0, r.nrows.pred:
        forO j, 0, r.ncols.pred:
          if i == j:
            op(r[i,j], x[i])
          else:
            op(r[i,j], 0)
  template `op MM`*(rr:untyped; xx:untyped):untyped =
    mixin op
    optimizeAst:
      subst(r,rr,x,xx,i,_,j,_):
        assert(r.nrows == x.nrows)
        assert(r.ncols == x.ncols)
        forO i, 0, r.nrows.pred:
          forO j, 0, r.ncols.pred:
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
      forO i, 0, r.len.pred:
        op(r[i], x[i], ty)
  template `op VSV`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,tx,_,i,_):
      assert(r.len == y.len)
      load(tx, x)
      forO i, 0, r.len.pred:
        op(r[i], tx, y[i])
  template `op VVV`*(r: typed; x,y: typed): untyped =
    mixin op
    #subst(r,rr,x,xx,y,yy,i,_):
    assert(r.len == y.len)
    assert(r.len == x.len)
    forO i, 0, r.len.pred:
      op(r[i], x[i], y[i])
  template `op MSS`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,tx,_,ty,_,i,_,j,_):
      assert(r.nrows == r.ncols)
      load(tx, x)
      load(ty, y)
      forO i, 0, r.nrows.pred:
        forO j, 0, r.ncols.pred:
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
      forO i, 0, r.nrows.pred:
        forO j, 0, r.ncols.pred:
          if i == j:
            op(r[i,j], x[i], ty)
          else:
            op(r[i,j], 0, 0)
  template `op MSV`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,tx,_,i,_,j,_):
      assert(r.nrows == y.len)
      assert(r.ncols == y.len)
      load(tx, x)
      forO i, 0, r.nrows.pred:
        forO j, 0, r.ncols.pred:
          if i == j:
            op(r[i,j], tx, y[i])
          else:
            op(r[i,j], 0, 0)
  template `op MVV`*(rr:typed; xx,yy:typed):untyped =
    subst(r,rr,x,xx,y,yy,i,_,j,_):
      assert(r.nrows == x.len)
      assert(r.ncols == x.len)
      assert(x.len == y.len)
      forO i, 0, r.nrows.pred:
        forO j, 0, r.ncols.pred:
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
      forO i, 0, r.nrows.pred:
        forO j, 0, r.ncols.pred:
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
      forO i, 0, r.nrows.pred:
        forO j, 0, r.ncols.pred:
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
      forO i, 0, r.nrows.pred:
        forO j, 0, r.ncols.pred:
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
      forO i, 0, r.nrows.pred:
        forO j, 0, r.ncols.pred:
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
      forO i, 0, r.nrows.pred:
        forO j, 0, r.ncols.pred:
          op(r[i,j], x[i,j], y[i,j])

makeMap2(add)
makeMap2(sub)

template imulVS*(rr:untyped; xx:untyped):untyped =
  mixin imul
  subst(r,rr,x,xx,i,_,tx,_):
    load(tx, x)
    forO i, 0, r.len.pred:
      imul(r[i], tx)

template imulMS*(r: typed; x: typed): untyped =
  mixin imul
  #subst(r,rr,x,xx,i,_,j,_,tx,_):
  #load(tx, x)
  let tx = x
  forO i, 0, r.nrows.pred:
    forO j, 0, r.ncols.pred:
      imul(r[i,j], tx)

template mulSVV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,i,_):
    mixin mul, imadd, assign
    assert(x.len == y.len)
    tmpvar(tr, r)
    mul(tr, x[0], y[0])
    forO i, 1, x.len.pred:
      imadd(tr, x[i], y[i])
    assign(r, tr)

template mulVVS*(r:typed; x,y:typed):untyped =
  mixin mul
  assert(r.len == x.len)
  load(ty, y)
  forO i, 0, r.len.pred:
    mul(r[i], x[i], ty)

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
]#
template mulVSVU*(r: typed; x,y: typed): untyped =
  mixin mul, `:=`
  assert(r.len == y.len)
  forO i, 0, r.len.pred:
    #mul(r[i], x, y[i])
    r[i] := mul(x, y[i])
template mulVSV*(r: typed; x,y: typed): untyped =
  # prepMultipleAccess(x)
  # prepLinearAccess(y)
  flattenCallArgs(mulVSVU, r, x, y)

template mulMMS*(rr:untyped; xx,yy:untyped):untyped =
  subst(r,rr,x,xx,y,yy,tx,_,i,_,j,_):
    #mixin mul
    assert(r.nrows == x.nrows)
    assert(r.ncols == x.ncols)
    load(ty, y)
    forO i, 0, r.nrows.pred:
      forO j, 0, r.ncols.pred:
        mul(r[i,j], x[i,j], ty)

template mulMSM*(rr:typed; xx,yy:typed):untyped =
  mixin mul
  subst(r,rr,x,xx,y,yy,tx,_,i,_,j,_):
    assert(r.nrows == y.nrows)
    assert(r.ncols == y.ncols)
    load(tx, x)
    forO i, 0, r.nrows.pred:
      forO j, 0, r.ncols.pred:
        mul(r[i,j], tx, y[i,j])

#[
template mulVMV*(rr: typed; xx,yy: typed): untyped =
  #[
  subst(r,rr,x,xx,y,yy,tr,_,ty,_,i,_,j,_,ty0r,_,ty0i,_,tyjr,_,tyji,_):
    mixin nrows, ncols, mul, imadd, assign, load1
    assert(x.nrows == r.len)
    assert(x.ncols == y.len)
    #when true:
    when false:
      tmpvar(tr, r)
      #var tr{.noInit.}:array[r.len,type(load1(r[0]))]
      load(ty0r, y[0].re)
      forO i, 0, x.nrows.pred:
        mulCCR(tr[i], x[i,0], ty0r)
      load(ty0i, y[0].im)
      forO i, 0, x.nrows.pred:
        imaddCCI(tr[i], x[i,0], ty0i)
      forO j, 1, x.ncols.pred:
        load(tyjr, y[j].re)
        forO i, 0, x.nrows.pred:
          imaddCCR(tr[i], x[i,j], tyjr)
        load(tyji, y[j].im)
        forO i, 0, x.nrows.pred:
          imaddCCI(tr[i], x[i,j], tyji)
      assign(r, tr)
      #forO i, 0, r.len.pred: assign(r[i], tr[i])
    else:
      tmpvar(tr, r)
      block:
        load(ty, y[0])
        forO i, 0, x.nrows.pred:
          mul(tr[i], x[i,0], ty)
      forO j, 1, x.ncols.pred:
        load(ty, y[j])
        forO i, 0, x.nrows.pred:
          imadd(tr[i], x[i,j], ty)
      assign(r, tr)
  ]#
  mixin nrows, ncols, mul, imadd, assign, load1
  template r: untyped = rr
  let x = xx
  let y = yy
  assert(x.nrows == r.len)
  assert(x.ncols == y.len)
  #tmpvar(tr_mulVMV, r)
  var tr_mulVMV{.noInit.}: type(load1(r))
  block:
    #let ty_mulVMV = y[0]
    let ty_mulVMV = load1(y[0])
    forO i, 0, x.nrows.pred:
      mul(tr_mulVMV[i], x[i,0], ty_mulVMV)
  forO j, 1, x.ncols.pred:
    #let ty_mulVMV = y[j]
    let ty_mulVMV = load1(y[j])
    forO i, 0, x.nrows.pred:
      imadd(tr_mulVMV[i], x[i,j], ty_mulVMV)
  assign(r, tr_mulVMV)
]#
#[
template mulVMVU*(r: typed; x,y: typed): untyped =
  mixin nrows, ncols, mul, imadd, assign, load1
  assert(x.nrows == r.len)
  assert(x.ncols == y.len)
  block:
    #let ty_mulVMV = y[0]
    let ty_mulVMV = load1(y[0])
    forO i, 0, x.nrows.pred:
      mul(r[i], x[i,0], ty_mulVMV)
  forO j, 1, x.ncols.pred:
    #let ty_mulVMV = y[j]
    let ty_mulVMV = load1(y[j])
    forO i, 0, x.nrows.pred:
      imadd(r[i], x[i,j], ty_mulVMV)
]#
template mulVMVU*(r: typed; x,y: typed): untyped =
  mixin nrows, ncols, mul, imadd, assign, `:=`
  assert(x.nrows == r.len)
  assert(x.ncols == y.len)
  forO i, 0, x.nrows.pred:
    var t_mulVMVU{.noInit.}: type(x[i,0]*y[0])
    t_mulVMVU := x[i,0] * y[0]
    forO j, 1, x.ncols.pred:
      imadd(t_mulVMVU, x[i,j], y[j])
    r[i] := t_mulVMVU
template mulVMV*(r: typed; x,y: typed): untyped =
  flattenCallArgs(mulVMVU, r, x, y)

template mulMMM*(rr:typed; xx,yy:typed):untyped =
  #[
  subst(r,rr,x,xx,y,yy,tr,_,i,_,j,_,k,_,txi0r,_,txi0i,_,txikr,_,txiki,_):
    mixin nrows, ncols, mul, imadd, assign, load1
    assert(r.nrows == x.nrows)
    assert(r.ncols == y.ncols)
    assert(x.ncols == y.nrows)
    var tr{.noInit.}:VectorArray[r.ncols,type(x[0,0]*y[0,0])]
    forO i, 0, r.nrows.pred:
      load(txi0r, x[i,0].re)
      forO j, 0, r.ncols.pred:
        mulCRC(tr[j], txi0r, y[0,j])
      load(txi0i, x[i,0].im)
      forO j, 0, r.ncols.pred:
        imaddCIC(tr[j], txi0i, y[0,j])
      forO k, 1, x.ncols.pred:
        load(txikr, x[i,k].re)
        forO j, 0, r.ncols.pred:
          imaddCRC(tr[j], txikr, y[k,j])
        load(txiki, x[i,k].im)
        forO j, 0, r.ncols.pred:
          imaddCIC(tr[j], txiki, y[k,j])
      forO j, 0, r.ncols.pred:
        assign(r[i,j], tr[j])
  ]#
  XoptimizeAst:
    subst(r,rr,x,xx,y,yy,tr,_,i,_,j,_,k,_,txi0r,_,txi0i,_,txikr,_,txiki,_):
      assert(x.nrows == r.nrows)
      assert(x.ncols == y.nrows)
      assert(r.ncols == y.ncols)
      mixin mul, imadd
      forO i, 0, r.nrows.pred:
        var tr{.noInit.}:VectorArray[r.ncols,type(x[0,0]*y[0,0])]
        #load(txi0, x[i,0])
        let txi0 = x[i,0]
        forO j, 0, r.ncols.pred:
          mul(tr[j], txi0, y[0,j])
        forO k, 1, x.ncols.pred:
          #load(txik, x[i,k])
          let txik = x[i,k]
          forO j, 0, r.ncols.pred:
            imadd(tr[j], txik, y[k,j])
        forO j, 0, r.ncols.pred:
          assign(r[i,j], tr[j])


template imaddSVV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,i,_):
    mixin imadd, assign
    assert(x.len == y.len)
    load(tr, r)
    forO i, 0, x.len.pred:
      imadd(tr, x[i], y[i])
    assign(r, tr)

#[
template imaddVSV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,i,_):
    mixin imadd, assign
    assert(r.len == y.len)
    load(tr, r)
    forO i, 0, r.len.pred:
      imadd(tr[i], x, y[i])
    assign(r, tr)
]#
template imaddVSV*(r: typed; xx,yy: typed): untyped =
  mixin imadd, assign
  let x = xx
  let y = yy
  assert(r.len == y.len)
  load(tr, r)
  forO i, 0, r.len.pred:
    imadd(tr[i], x, y[i])
  assign(r, tr)

#[
template imaddVMV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,ty,_,i,_,j,_,tyjr,_,tyji,_):
    mixin nrows, ncols, mul, imadd, assign, load1
    assert(x.nrows == r.len)
    assert(x.ncols == y.len)
    when true:
    #when false:
      load(tr, r)
      #var tr{.noInit.}:array[r.len,type(load1(r[0]))]
      #forO i, 0, r.len.pred: assign(tr[i], r[i])
      forO j, 0, x.ncols.pred:
        load(tyjr, y[j].re)
        forO i, 0, x.nrows.pred:
          imaddCCR(tr[i], x[i,j], tyjr)
        load(tyji, y[j].im)
        forO i, 0, x.nrows.pred:
          imaddCCI(tr[i], x[i,j], tyji)
      assign(r, tr)
      #forO i, 0, r.len.pred: assign(r[i], tr[i])
    else:
      load(tr, r)
      forO j, 0, x.ncols.pred:
        load(tyr, asReal(y[j].re))
        forO i, 0, x.nrows.pred:
          imadd(tr[i], x[i,j], tyr)
        load(tyi, asImag(y[j].im))
        forO i, 0, x.nrows.pred:
          imadd(tr[i], x[i,j], tyi)
      assign(r, tr)
]#
template imaddVMVU*(r: typed; x,y: typed): untyped =
  mixin nrows, ncols, imadd, `:=`
  assert(x.nrows == r.len)
  assert(x.ncols == y.len)
  forO i, 0, x.nrows.pred:
    #var t = x[i,0] * y[0]
    #forO j, 1, x.ncols.pred:
    #  imadd(t, x[i,j], y[j])
    #r[i] += t
    var t = r[i]
    forO j, 0, x.ncols.pred:
      imadd(t, x[i,j], y[j])
    r[i] := t
template imaddVMV*(r: typed; x,y: typed): untyped =
  flattenCallArgs(imaddVMVU, r, x, y)

template imaddMMM*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,ty,_,i,_,j,_,k,_,tyjr,_,tyji,_):
    mixin nrows, ncols, mul, imadd, assign, load1
    assert(r.nrows == x.nrows)
    assert(r.ncols == y.ncols)
    assert(x.ncols == y.nrows)
    #when true:
    when false:
      load(tr, r)
      forO i, 0, r.nrows.pred:
        forO k, 0, x.ncols.pred:
          load(txikr, x[i,k].re)
          forO j, 0, r.ncols.pred:
            imaddCRC(tr[i,j], txikr, y[k,j])
          load(txiki, x[i,k].im)
          forO j, 0, r.ncols.pred:
            imaddCIC(tr[i,j], txiki, y[k,j])
      assign(r, tr)
    else:
      forO i, 0, r.nrows.pred:
        var tr{.noInit.}:VectorArray[r.ncols,type(x[0,0]*y[0,0])]
        forO j, 0, r.ncols.pred:
          assign(tr[j], r[i,j])
        forO k, 0, x.ncols.pred:
          load(txikr, x[i,k].re)
          forO j, 0, r.ncols.pred:
            imaddCRC(tr[j], txikr, y[k,j])
          load(txiki, x[i,k].im)
          forO j, 0, r.ncols.pred:
            imaddCIC(tr[j], txiki, y[k,j])
        forO j, 0, r.ncols.pred:
          assign(r[i,j], tr[j])

template imsubVSV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tx,_,i,_):
    mixin imsub
    assert(r.len == y.len)
    load(tx, x)
    forO i, 0, r.len.pred:
      imsub(r[i], x, y[i])

template imsubVMV*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,ty,_,i,_,j,_):
    mixin imsub
    assert(x.nrows == r.len)
    assert(x.ncols == y.len)
    when true:
      load(tr, r)
      forO j, 0, x.ncols.pred:
        load(ty, y[j])
        forO i, 0, x.nrows.pred:
          imsub(tr[i], x[i,j], ty)
      assign(r, tr)
    else:
      load(tr, r)
      forO j, 0, x.ncols.pred:
        load(tyr, asReal(y[j].re))
        forO i, 0, x.nrows.pred:
          imsub(tr[i], x[i,j], tyr)
        load(tyi, asImag(y[j].im))
        forO i, 0, x.nrows.pred:
          imsub(tr[i], x[i,j], tyi)
      assign(r, tr)

template imsubMMM*(rr:typed; xx,yy:typed):untyped =
  subst(r,rr,x,xx,y,yy,tr,_,ty,_,i,_,j,_,k,_,txikr,_,txiki,_):
    mixin nrows, ncols, imsubCRC, imsubCIC, assign, load1
    assert(r.nrows == x.nrows)
    assert(r.ncols == y.ncols)
    assert(x.ncols == y.nrows)
    load(tr, r)
    forO i, 0, r.nrows.pred:
      forO k, 0, x.ncols.pred:
        load(txikr, x[i,k].re)
        forO j, 0, r.ncols.pred:
          imsubCRC(tr[i,j], txikr, y[k,j])
        load(txiki, x[i,k].im)
        forO j, 0, r.ncols.pred:
          imsubCIC(tr[i,j], txiki, y[k,j])
    assign(r, tr)

template msubVSVV*(rr:typed; xx,yy,zz:typed):untyped =
  subst(r,rr,x,xx,y,yy,z,zz,i,_):
    mixin msub
    assert(r.len == y.len)
    assert(r.len == z.len)
    forO i, 0, r.len.pred:
      msub(r[i], x, y[i], z[i])
