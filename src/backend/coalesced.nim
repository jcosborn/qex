#[

Following Nvidia's idea of coalesced_ptr in C++, we use a wrapper
object type to hide the actual coalesced memory layout here.
Original comments from Nvidia's coalesced_ptr.h follows:

  A smart pointer that automatically provide coalesced memory
  transcations for arrays of arbtrary structures.  Given a structure
  T, of size S bytes, e.g.,

  struct T {
    char a[S];
  }

  in an array with sites elements

  T t[sites];

  using a coalesced_ptr will split the structure for reading and
  writing to memory as an array of structures of array of structures (AoSoAoS),
  where:
    - the inner structure size is given by memory_word_size
    - the inner array size is given by site_vector
    - the outer structure size is given by sizeof(T)/memory_word_size
    - the outer array size is given by sites/site_vector

]#

import macros, base/metaUtils

type
  Coalesced*[V,M:static[int],T] = object
    ## `V`: Inner array length.
    ## `M`: Number of RegisterWords in a MemoryWord, the granularity of memory transactions.
    p*: ptr T                   ## pointer to an array of T
    n*: int                     ## the length of the array being coalesced
  CoalescedObj[V,M:static[int],T] = object
    o*: Coalesced[V,M,T]
    i*: int                     # the index to which we asks

const llbits = currentSourcePath()[0..^14] & "llbits.h"
type
  RegisterWord* {.importc, header:llbits.} = object # Word fits in a register, 4 bytes for current GPU
  MemoryWord1 {.importc, header:llbits.} = object # Word fits in a register, 4 bytes for current GPU
    a*: array[1,RegisterWord]
  MemoryWord2 {.importc, header:llbits.} = object # Word fits in a register, 4 bytes for current GPU
    a*: array[2,RegisterWord]
  MemoryWord4 {.importc, header:llbits.} = object # Word fits in a register, 4 bytes for current GPU
    a*: array[4,RegisterWord]
  MemoryWord8 {.importc, header:llbits.} = object # Word fits in a register, 4 bytes for current GPU
    a*: array[8,RegisterWord]
  MemoryWord16 {.importc, header:llbits.} = object # Word fits in a register, 4 bytes for current GPU
    a*: array[16,RegisterWord]
  MemoryWord32 {.importc, header:llbits.} = object # Word fits in a register, 4 bytes for current GPU
    a*: array[32,RegisterWord]
template MemoryWord(M:static[int]):untyped =
  when 1 == M: MemoryWord1
  elif 2 == M: MemoryWord2
  elif 4 == M: MemoryWord4
  elif 8 == M: MemoryWord8
  elif 16 == M: MemoryWord16
  elif 32 == M: MemoryWord32

template sizeOf*(t:typedesc[RegisterWord]):int = 4

# Nim doesn't know the size of any struct for sure without the help of a C/C++ compiler.
# So we use a C++ compiler to check if the user has provided a correct size.
# The following C++ code only works with c++11 or later.
{.emit:"""
#if __cplusplus >= 201103L
  template <typename ToCheck, std::size_t ProvidedSize, std::size_t RealSize = sizeof(ToCheck)>
  void coalesced_check_size() {static_assert(ProvidedSize == RealSize, "newCoalesced got the wrong size!");}
#else
  #define coalesced_check_size(type,size) typedef char ProvidedWrongSizeForType##type[2*!!(sizeof(type)==(size))-1]
#endif
""".}
template getSize*(T:typedesc):untyped =
  when compiles((const size = sizeof(T))):
    const size = sizeof(T)
  else:
    mixin structSize
    const size = structSize(T)
  {.emit:"""
    #if __cplusplus >= 201103L
      coalesced_check_size<`T`,`size`>();
    #else
      coalesced_check_size(`T`,`size`);
    #endif
  """.}
  size

# Nim bug as of 8/7/2017, cannot overload init/newCoalesced.
# Overloaded type matching would SIGSEGV.
proc initCoalesced*[V,M:static[int],T](x:var Coalesced[V,M,T], p:ptr T, n:int) =
  const
    size = getSize(T)
    N = size div (M*sizeof(RegisterWord))
  when N*(M*sizeof(RegisterWord)) != size: {.fatal:"sizeof(T) must be divisible by memory word size."}
  if n mod V != 0:
    echo "Array length for Coalesced must be multiples of V = ",V
    quit 1
  x.p = p
  x.n = n
proc newCoalesced*[T](V,M:static[int], p:ptr T, n:int):auto {.noinit.} =
  var r {.noinit.}:Coalesced[V,M,T]
  r.initCoalesced(p,n)
  r

template `[]`*(x:Coalesced, ix:int):untyped = CoalescedObj[x.V,x.M,x.T](o:x, i:ix)
template len*(x:Coalesced):untyped = x.n

template fromCoalesced*(x:CoalescedObj):untyped =
  const N = getSize(x.T) div (x.M*sizeof(RegisterWord))
  type A = ptr Uncheckedarray[MemoryWord(x.M)]
  var r {.noinit.}: x.T
  let offset = (x.i div x.V)*N*x.V + x.i mod x.V
  unrollfor:
    for j in 0..N-1: cast[A](r.addr)[j] = cast[A](x.o.p)[offset + j*x.V]
  r
macro `[]`*(x:CoalescedObj, ys:varargs[untyped]):untyped =
  let o = newCall(bindsym"fromCoalesced", x)
  if ys.len == 0:
    result = o
  else:
    result = newCall("[]", o)
    for y in ys: result.add y

proc `:=`*[V,M:static[int],X,Y](x:CoalescedObj[V,M,X], y:var Y) {.inline.} =
  when Y is x.T:
    const N = getSize(x.T) div (x.M*sizeof(RegisterWord))
    type A = ptr Uncheckedarray[MemoryWord(x.M)]
    let offset = (x.i div x.V)*N*x.V + x.i mod x.V
    unrollfor:
      for j in 0..N-1: cast[A](x.o.p)[offset + j*x.V] = cast[A](y.addr)[j]
  else:
    mixin `:=`
    var ty {.noinit.}:x.T
    ty := y
    x := ty
template `:=`*[V,M:static[int],X,Y](x:CoalescedObj[V,M,X], y:Y) =
  mixin `:=`
  var ty {.noinit.}:x.T
  ty := y
  x := ty

proc `*`*[VX,MX,VY,MY:static[int],X,Y](x:CoalescedObj[VX,MX,X], y:CoalescedObj[VY,MY,Y]):auto {.noinit,inline.} =
  let
    tx {.noinit.} = fromCoalesced(x)
    ty {.noinit.} = fromCoalesced(y)
  mixin `*`
  tx * ty

template `+=`*[Y](x:CoalescedObj, y:Y) = x := fromCoalesced(x) + y

when isMainModule:
  import strutils
  type T = array[6,int32]
  proc structSize(t:typedesc[T]):int = 24
  var x {.noinit.}: array[16,T]
  let p = newCoalesced(8, 2, x[0].addr, x.len)
  # var p:Coalesced[8,3,T]
  # p.initCoalesced(x[0].addr, x.len)
  for i in 0..<p.len:
    var t {.noinit.}: T
    for j in 0..<t.len: t[j] = int32(100*i + j)
    p[i] := t
  var s:string
  s = "Lexical order: p = {"
  for i in 0..<p.len:
    let t = p[i][]
    s &= "\n["
    for j in 0..<t.len: s &= " " & align($t[j],4)
    s &= " ]"
  s &= "}"
  echo s
  s = "Memory layout: x = {"
  var c = 0
  for i in 0..<x.len:
    for j in 0..<x[0].len:
      if 0 == c mod (p.V*p.M): s &= "\n"
      inc c
      s &= " " & align($x[i][j],4)
  s &= "}"
  echo s
