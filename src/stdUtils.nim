import macros
import strUtils
import metaUtils

type
  cArray*{.unchecked.}[T] = array[0,T]

template ptrInt*(x:untyped):untyped = cast[ByteAddress](x)
template addrInt*(x:untyped):untyped = cast[ByteAddress](addr(x))
template unsafeAddrInt*(x:untyped):untyped = cast[ByteAddress](addr(x))

template `$&`*(x: untyped): string =
  toHex(unsafeAddrInt(x))

proc `|`*(x: int, d: int): string =
  result = $x
  var plen = d.abs-len(result)
  if plen<0: plen = 0
  let pad = repeat(' ', plen)
  if d >= 0:
    result = pad & result
  else:
    result = result & pad
proc `|`*(s: string, d: int): string =
  let pad = repeat(' ', d.abs-len(s))
  if d >= 0:
    result = pad & s
  else:
    result = s & pad
proc `|`*(s: string, d: tuple[w:int,c:char]): string =
  let pad = repeat(d.c, d.w.abs-len(s))
  if d.w >= 0:
    result = pad & s
  else:
    result = s & pad
proc `|`*(f: float, d: tuple[w,p: int]): string =
  result = formatFloat(f, ffDecimal, d.p)
  let p = d.w.abs-len(result)
  let pad = if p>0: repeat(' ', p) else: ""
  if d.w >= 0:
    result = pad & result
  else:
    result = result & pad
proc `|`*(f: float, d: int): string =
  $f | d

proc `*`*[T](x:openArray[T], y:int):auto {.inline.} =
  let n = x.len
  var r:array[n,T]
  for i in 0..<n:
    r[i] = x[i]
  r
proc `+`*[T](x,y:openArray[T]):auto {.inline.} =
  let n = x.len
  var r:array[n,T]
  for i in 0..<n:
    r[i] = x[i] + y[i]
  r
proc `+=`*[T](x:var openArray[T], y: openArray[T]) {.inline.} =
  let n = x.len
  for i in 0..<n:
    x[i] += y[i]
template makeArrayOverloads(n:int):untyped =
  proc `+`*[T](x,y:array[n,T]):array[n,T] {.inline.} =
    for i in 0..<x.len:
      result[i] = x[i] + y[i]
  proc `*`*[T](x:array[n,T], y:int):array[n,T] {.inline.} =
    for i in 0..<x.len:
      result[i] = x[i] * T(y)
  proc `:=`*[T1,T2](r:var array[n,T1]; x:array[n,T2]) =
    for i in 0..<r.len:
      r[i] = T1(x[i])
makeArrayOverloads(4)
makeArrayOverloads(8)
makeArrayOverloads(16)

proc sum*[T](x:openArray[T]): auto =
  result = x[0]
  for i in 1..<x.len: result += x[i]

macro echoImm*(s:varargs[typed]):auto =
  result = newEmptyNode()
  #echo s.treeRepr
  var t = ""
  for c in s.children():
    if c.kind == nnkStrLit:
      t &= c.strVal
    else:
      t &= c.toStrLit.strVal
  echo t

template ctrace* =
  const ii = instantiationInfo()
  echoImm "ctrace: ", ii

template declareVla(v,t,n:untyped):untyped =
  type Vla{.gensym.} = distinct t
  #var v{.noInit,codeGenDecl:"$# $#[" & n.astToStr & "]".}:Vla
  #var v{.noInit,codeGenDecl:"$# $#[`n`]".}:Vla
  var v{.noInit,noDecl.}:Vla
  {.emit:"`Vla` `v`[`n`];".}
  template len(x:Vla):untyped = n
  template `[]`(x:Vla; i:untyped):untyped =
    (cast[ptr cArray[t]](unsafeAddr(x)))[][i]
  template `[]=`(x:var Vla; i,y:untyped):untyped =
    (cast[ptr cArray[t]](addr(x)))[][i] = y

proc `$`*[T](x:openArray[T]):string =
  var t = newSeq[string]()
  var len = 0
  for e in x:
    let s = $e
    t.add(s)
    len += s.len
  #echo len
  #echo t[0]
  if len < 60:
    result = t.join(" ")
  else:
    result = ""
    for i,v in t:
      result &= ($i & ":" & v & "\n")

macro toLit*(s:static[string]):auto =
  result = newLit(s)

template warn*(s:varargs[string,`$`]) =
  let ii = instantiationInfo()
  echo "warning (", ii.filename, ":", ii.line, "):"
  echo "  ", s.join
  
#[
when isMainModule:
  proc test(n:int) =
    declareVla(x, float, n)
    let n2 = n div 2
    block:
      declareVla(y, float, n2)
      #{.emit:"""printf("%p\n", &x[0]);""".}
      x[0] = 1
      echo x[0]
      echo x.len
      echo y.len
  test(10)
  test(20)
]#
