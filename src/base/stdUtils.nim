import macros
import strUtils
import metaUtils
import os

type
  cArray*[T] = UncheckedArray[T]
template `[]`*(x: cArray): untyped = addr x[0]
template `&`*(x: ptr cArray): untyped = addr x[0]

template ptrInt*(x:untyped):untyped = cast[ByteAddress](x)
template addrInt*(x:untyped):untyped = cast[ByteAddress](addr(x))
template unsafeAddrInt*(x:untyped):untyped = cast[ByteAddress](addr(x))

proc newSeqU*[T](n: int): seq[T] =
  result = newSeqOfCap[T](n)
  result.setLen(n)

iterator range*[T: SomeInteger](count: T): T =
  var res = T(0)
  while res < count:
    yield res
    inc res

var paramNames = newSeq[string](0)
var paramValues = newSeq[string](0)

proc addParam(s,r: string) =
  var i = paramNames.find(s)
  if i>=0:
    paramValues[i] = r
  else:
    paramNames.add s
    paramValues.add r

template echoParams*() =
  mixin echo
  for i in 0..<paramNames.len:
    echo paramNames[i], ": ", paramValues[i]

template makeTypeParam(name,typ,deflt,cnvrt: untyped): untyped {.dirty.} =
  proc name*(s: string, d=deflt): typ =
    result = d
    let n = paramCount()
    for i in 1..n:
      let p = paramstr(i)
      if p.startsWith('-'&s&':'):
        let ll = s.len + 2
        result = cnvrt(p[ll..^1])
    addParam(s, $result)

makeTypeParam(intParam, int, 0, parseInt)
makeTypeParam(floatParam, float, 0.0, parseFloat)
makeTypeParam(strParam, string, "", string)
template stringParam*(x,y: untyped): untyped = strParam(x,y)

proc intSeqParam*(s: string, d: seq[int] = @[]): seq[int] =
  result = d
  let n = paramCount()
  for i in 1..n:
    let p = paramstr(i)
    if p.startsWith('-'&s&':'):
      result.setLen(0)
      let ll = s.len + 2
      for c in split(p[ll..^1], ','):
        result.add parseInt(c)
  addParam(s, join(result," "))

proc floatSeqParam*(s: string, d: seq[float] = @[]): seq[float] =
  result = d
  let n = paramCount()
  for i in 1..n:
    let p = paramstr(i)
    if p.startsWith('-'&s&':'):
      result.setLen(0)
      let ll = s.len + 2
      for c in split(p[ll..^1], ','):
        result.add parseFloat(c)
  addParam(s, join(result," "))

template CLIset*(p:typed, n:untyped, prefix:string, runifset:untyped) =
  mixin echo
  let
    o = p.n
    s = prefix & astToStr(n)
  when compiles(strParam(s, p.n)):
    p.n = type(p.n)strParam(s, p.n)
  elif compiles(intParam(s, p.n)):
    p.n = type(p.n)intParam(s, p.n)
  elif compiles(floatParam(s, p.n)):
    p.n = type(p.n)floatParam(s, p.n)
  elif compiles(intSeqParam(s, p.n)):
    p.n = type(p.n)intSeqParam(s, p.n)
  elif compiles(floatSeqParam(s, p.n)):
    p.n = type(p.n)floatSeqParam(s, p.n)
  else:
    {.fatal:"Cannot set argument "&s&" of "&astToStr(p)&" for command line.".}
  if o != p.n:
    runifset
    echo "Customize $# : $# -> $#"%[s, $o, $p.n]
template CLIset*(p:typed, n:untyped, prefix = "") =
  p.CLIset n, prefix:
    discard

template `$&`*(x: untyped): string =
  toHex(unsafeAddrInt(x))

proc `|`*(s: string, d: tuple[w:int,c:char]): string =
  let p = abs(d.w) - len(s)
  let pad = if p>0: repeat(d.c, p) else: ""
  if d.w >= 0:
    result = pad & s
  else:
    result = s & pad
proc `|`*(s: string, d: int): string =
  s | (d,' ')
proc `|`*(x: int, d: int): string =
  ($x) | d
proc `|`*(f: float, d: tuple[w,p: int]): string =
  if d.p<0:
    formatFloat(f, ffDecimal, -d.p) | d.w
  else:
    formatFloat(f, ffDefault, d.p) | d.w
proc `|`*(f: float, d: int): string =
  f | (d,d)
template `|-`*(x:SomeNumber, y: int): untyped =
  x | -y

proc indexOf*[T](x: openArray[T], y: any): int =
  let n = x.len
  while result<n and x[result]!=y: inc result

proc `+`*(x: SomeNumber, y: array): auto {.inline,noInit.} =
  var r: array[y.len, type(x+y[0])]
  for i in 0..<r.len:
    r[i] = x + y[i]
  r

proc `*`*(x: array, y: SomeNumber): auto {.inline,noInit.} =
  var r: array[x.len, type(x[0]*y)]
  for i in 0..<r.len:
    r[i] = x[i] * y
  r
proc `*`*(x: SomeNumber, y: array): auto {.inline,noInit.} =
  var r: array[y.len, type(x*y[0])]
  for i in 0..<r.len:
    r[i] = x * y[i]
  r
proc `*`*[T](x: SomeNumber, y: seq[T]): seq[T] {.inline,noInit.} =
  result.newSeq(y.len)
  for i in 0..<result.len:
    result[i] = x * y[i]

proc `/`*(x: SomeNumber, y: array): auto {.inline,noInit.} =
  var r: array[y.len, type(x/y[0])]
  for i in 0..<r.len:
    r[i] = x / y[i]
  r

#proc `+`*[A:array](x,y: A): A {.inline,noInit.} =
#  for i in 0..<result.len:
#    result[i] = x[i] + y[i]

proc `:=`*[N,T1,T2](r: var array[N,T1], x: array[N,T2]) {.inline.} =
  mixin `:=`
  const n = r.len
  for i in 0..<n:
    r[i] := x[i]

proc `+`*[N,T1,T2](x: array[N,T1], y: array[N,T2]): auto {.inline,noInit.} =
  const n = x.len
  var r: array[n, type(x[0]+y[0])]
  for i in 0..<n:
    r[i] = x[i] + y[i]
  r

proc `-`*[N,T1,T2](x: array[N,T1], y: array[N,T2]): auto {.inline,noInit.} =
  const n = x.len
  var r: array[n, type(x[0]-y[0])]
  for i in 0..<n:
    r[i] = x[i] - y[i]
  r

proc `*`*[N,T1,T2](x: array[N,T1], y: array[N,T2]): auto {.inline,noInit.} =
  const n = x.len
  var r: array[n, type(x[0]*y[0])]
  for i in 0..<n:
    r[i] = x[i] * y[i]
  r

proc `+=`*[T](r: var openArray[T], x: openArray[T]) {.inline.} =
  let n = r.len
  for i in 0..<n:
    r[i] += x[i]

#[
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
]#

#proc sum*[T](x: openArray[T]): T =
#  for i in 0..<x.len: result += x[i]

proc product*[T](x: openArray[T]): T =
  result = T(1)
  for i in 0..<x.len: result *= x[i]

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

#[
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
]#

macro toLit*(s:static[string]):auto =
  result = newLit(s)

template warn*(s:varargs[string,`$`]) =
  let ii = instantiationInfo()
  echo "warning (", ii.filename, ":", ii.line, "):"
  echo "  ", s.join

proc factor*(n: int): seq[int] =
  result.newSeq(0)
  var x = n
  if x<0:
    result.add(-1)
    x = -x
  if x<2: result.add x
  while x>1:
    var k = 2
    if (x and 1) != 0:
      k = 3
      while (x mod k) != 0: k += 2
    result.add k
    x = x div k


when isMainModule:
  #[
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

  template testFactor(n: int) =
    echo "factor(", n, ") = ", factor(n)
  for i in -2..20:
    testFactor(i)
