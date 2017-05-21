import strUtils
import stdUtils

type
  alignedMem*[T] = object
    len*: int
    align*: int
    stride*: int
    bytes*: int
    mem*: ref cArray[char]
    data*: ptr cArray[T]

proc unsafeNewU*[T](a: var ref T, size: Natural) =
  {.emit: "N_NIMCALL(void*, newObjNoInit)(TNimType* typ0, NI size0);".}
  {.emit: "#define newObj newObjNoInit".}
  unsafeNew(a, size)
  {.emit: "#undef newObj".}

proc ptrAlign[T](p:ptr T; a:int):ptr T =
  let x = cast[ByteAddress](p)
  let a1 = a - 1
  let y = x + (a1-((x+a1) mod a))
  #echo x, ":", y
  result = cast[type(result)](y)

proc new*[T](t:var alignedMem[T], n:int, align:int=64) =
  t.len = n
  t.align = align
  t.stride = sizeof(T)
  t.bytes = t.len * t.stride + t.align
  unsafeNew(t.mem, t.bytes)
  t.data = ptrAlign(cast[ptr cArray[T]](t.mem[0].addr), align)
proc newU*[T](t:var alignedMem[T], n:int, align:int=64) =
  t.len = n
  t.align = align
  t.stride = sizeof(T)
  t.bytes = t.len * t.stride + t.align
  unsafeNewU(t.mem, t.bytes)
  t.data = ptrAlign(cast[ptr cArray[T]](t.mem[0].addr), align)
proc newAlignedMem*[T](t:var alignedMem[T], n:int, align:int=64) =
  new(t, n, align)
proc newAlignedMem*[T](n:int, align:int=64): alignedMem[T] =
  newAlignedMem[T](result, n, align)
proc newAlignedMemU*[T](t:var alignedMem[T], n:int, align:int=64) =
  newU(t, n, align)
proc newAlignedMemU*[T](n:int, align:int=64): alignedMem[T] =
  newAlignedMemU[T](result, n, align)

template low*(s:alignedMem):untyped = 0
template high*(s:alignedMem):untyped = s.len-1
proc `[]`*[T](s:alignedMem[T], i:SomeInteger):var T =
  result = s.data[i]
#template `[]`*[T](s:alignedMem[T], i:SomeInteger):untyped = s.data[i]
template `[]`*[T](s:var alignedMem[T], i:SomeInteger):untyped = s.data[i]
template `[]=`*[T](s:var alignedMem[T], i:SomeInteger, v:untyped) =
  s.data[i] = v

when isMainModule:
  var x: alignedMem[float]
  newAlignedMem(x, 10)
  let c0 = cast[ByteAddress](x.mem[0].addr)
  echo c0, " ", toHex(c0,8)
  let x0 = cast[ByteAddress](x[0].addr)
  echo x0, " ", toHex(x0,8)

  for i in x.low..x.high:
    x[i] = float(i)
  for i in 0..<x.len:
    assert(x[i] == float(i))
