import macros
import parseutils
import strUtils
import algorithm

type
  cArray*[T] = UncheckedArray[T]
#template `[]`*(x: cArray): untyped = addr x[0]
#template `&`*(x: ptr cArray): untyped = addr x[0]

template ptrInt*(x: auto): auto = cast[uint](x)
template addrInt*(x: auto): auto = cast[uint](addr(x))
template unsafeAddrInt*(x: auto): auto = cast[uint](addr(x))
template toHex*(x: ptr auto): auto = toHex(cast[uint](x))
template toHex*(x: pointer): auto = toHex(cast[uint](x))

type
  ConstInt* {.importc:"const int".} = object
proc constCast*(x: ptr ConstInt): ptr cint {.importc:"(int *)",nodecl.}

proc newSeqU*[T](n: int): seq[T] =
  result = newSeqOfCap[T](n)
  result.setLen(n)

iterator range*[T: SomeInteger](count: T): T =
  var res = T(0)
  while res < count:
    yield res
    inc res

proc isInteger*(s: string):bool =
  var t:int
  parseInt(s, t) == s.len

template `$&`*(x: auto): string =
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
proc `|`*(x: SomeInteger, d: int): string =
  ($x) | d
proc `|`*(f: SomeFloat, d: tuple[w,p: int]): string =
  if d.p<0:
    formatFloat(f, ffDecimal, -d.p) | d.w
  else:
    formatFloat(f, ffDefault, d.p) | d.w
proc `|`*(f: float, d: int): string =
  f | (d,d)
template `|-`*(x:SomeNumber, y: int): auto =
  x | -y

proc values*(e: typedesc[enum]): string =
  result = ""
  var sep = ""
  for v in e:
    result &= sep & $v
    sep = " "

proc indexOf*[T](x: openArray[T], y: auto): int =
  let n = x.len
  while result<n and x[result]!=y: inc result

proc `+`*[N,T](x: SomeNumber, y: array[N,T]): auto {.inline,noInit.} =
  var r: array[N, type(x+y[0])]
  for i in 0..<r.len:
    r[i] = x + y[i]
  r

proc `*`*[N,T](x: array[N,T], y: SomeNumber): auto {.inline,noInit.} =
  var r: array[N, type(x[0]*y)]
  for i in 0..<r.len:
    r[i] = x[i] * y
  r
proc `*`*[N,T](x: SomeNumber, y: array[N,T]): auto {.inline,noInit.} =
  var r: array[N, type(x*y[0])]
  for i in 0..<r.len:
    r[i] = x * y[i]
  r
proc `*`*[T](x: seq[T], y: SomeNumber): seq[T] {.inline,noInit.} =
  result.newSeq(x.len)
  for i in 0..<result.len:
    result[i] = x[i] * y
proc `*`*[T](x: SomeNumber, y: seq[T]): seq[T] {.inline,noInit.} =
  result.newSeq(y.len)
  for i in 0..<result.len:
    result[i] = x * y[i]

proc `/`*[N,T](x: SomeNumber, y: array[N,T]): auto {.inline,noInit.} =
  var r: array[N, type(x/y[0])]
  for i in 0..<r.len:
    r[i] = x / y[i]
  r

proc `mod`*[N,T](x: array[N,T], y: SomeNumber): auto {.inline,noInit.} =
  var r: array[N, type(x[0] mod y)]
  for i in 0..<r.len:
    r[i] = x[i] mod y
  r

#proc `+`*[A:array](x,y: A): A {.inline,noInit.} =
#  for i in 0..<result.len:
#    result[i] = x[i] + y[i]

proc assign*[N,T](r: var array[N,T], x: SomeNumber) {.inline.} =
  mixin assign
  const n = r.len
  for i in 0..<n:
    assign(r[i], x)

template `:=`*[N,T](r: var array[N,T], x: SomeNumber) =
  assign(r, x)

proc assign*[N1,T1,N2,T2](r: var array[N1,T1], x: array[N2,T2]) {.inline.} =
  mixin assign
  const n1 = r.len
  const n2 = x.len
  when n1 == n2:
    for i in 0..<n1:
      assign(r[i], x[i])
  elif n1 > n2:
    const nr = n1 div n2
    when nr*n2 != n1:
      static:
        echo "error: ':=' array sizes ", n1, " and ", n2
        quit()
    for i in 0..<n2:
      var ri = cast[ptr array[nr,T1]](addr r[nr*i])
      assign(ri[], x[i])

template `:=`*[N1,T1,N2,T2](r: var array[N1,T1], x: array[N2,T2]) =
  assign(r, x)

#proc `:=`*[N,T1,T2](r: var array[N,T1], x: array[N,T2]) {.inline.} =
#  mixin `:=`
#  const n = r.len
#  for i in 0..<n:
#    r[i] := x[i]
#template assign*[N,T1,T2](r: var array[N,T1], x: array[N,T2]) =
#  r := x

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

proc `+`*[T1,T2](x: seq[T1], y: openArray[T2]): auto {.inline,noInit.} =
  mixin evalType
  let n = min(x.len, y.len)
  var r = newSeq[evalType(x[0]+y[0])](n)
  for i in 0..<n:
    r[i] = x[i] + y[i]
  r

proc `-`*[T1,T2](x: seq[T1], y: openArray[T2]): auto {.inline,noInit.} =
  mixin evalType
  let n = min(x.len, y.len)
  var r = newSeq[evalType(x[0]+y[0])](n)
  for i in 0..<n:
    r[i] = x[i] - y[i]
  r

proc `:=`*[T](r: var openArray[T], x: openArray[T]) {.inline.} =
  let n = r.len
  for i in 0..<n:
    r[i] := x[i]

proc `+=`*[T](r: var openArray[T], x: openArray[T]) {.inline.} =
  let n = r.len
  for i in 0..<n:
    r[i] += x[i]

proc `-`*[T](x: seq[T]): seq[T] {.inline.} =
  let n = x.len
  result.newSeq(n)
  for i in 0..<n:
    result[i] = -x[i]

#proc sum*[T](x: openArray[T]): T =
#  for i in 0..<x.len: result += x[i]

proc product*[T](x: openArray[T]): T =
  result = T(1)
  for i in 0..<x.len: result *= x[i]

template memoize*(xi,yi,zi:int, computeStmt:untyped):auto =
  ## Simple nested seq based memoization for 3 int arguments.
  let
    x = xi
    y = yi
    z = zi
  type T = typeof(computeStmt)
  var cache {.global.} = newseq[seq[seq[T]]]()
  if (let lmax=cache.len; x>=lmax):
    cache.setlen(x+1)
    for l in lmax..x:
      cache[l] = newseq[seq[T]]()
  if (let lmax=cache[x].len; y>=lmax):
    cache[x].setlen(y+1)
    for l in lmax..y:
      cache[x][l] = newseq[T]()
  if z>=cache[x][y].len:
    cache[x][y].setlen(z+1)
  if cache[x][y][z]==nil:
    cache[x][y][z] = computeStmt
  cache[x][y][z]

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

#[
template declareVla(v,t,n:untyped):untyped =
  type Vla{.gensym.} = distinct t
  #var v{.noInit,codeGenDecl:"$# $# [" & n.astToStr & "]".}:Vla
  #var v{.noInit,codeGenDecl:"$# $# [`n`]".}:Vla
  var v{.noInit,noDecl.}:Vla
  {.emit:"`Vla` `v`[`n`];".}
  template len(x:Vla):untyped = n
  template `[]`(x:Vla; i:untyped):untyped =
    (cast[ptr cArray[t]](unsafeAddr(x)))[][i]
  template `[]=`(x:var Vla; i,y:untyped):untyped =
    (cast[ptr cArray[t]](addr(x)))[][i] = y
]#

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

proc getDivisors*[T](n: T): seq[T] =
  result.newSeq(1)
  result[0] = 1
  for i in 2..n:
    if n mod i == 0:
      result.add i

# scale local lattice (ll) by nrank to get global lattice (gl)
template latticeFromLocalLatticeImpl*(gl: var openArray, ll: openArray, nrank: int) =
  var fs = factor(nrank)
  reverse fs
  let nd = ll.len
  for i in 0..<nd: gl[i] = ll[i]
  var k = nd - 1
  for f in fs:
    gl[k] *= f
    k = (k+nd-1) mod nd

proc latticeFromLocalLattice*[T:seq|array](ll: T, nrank: int): T =
  when T is seq:
    result.newSeq(ll.len)
  latticeFromLocalLatticeImpl(result, ll, nrank)

macro rangeLow*(r: typedesc[range]):auto =
  #echo r.treerepr
  #echo r.getType.treerepr
  #echo r.getTypeImpl.treerepr
  #echo r.getTypeInst.treerepr
  #error("rangeLow(typedesc)")
  #return r
  let t = r.getTypeImpl
  t.expectKind(nnkBracketExpr)
  let t1 = t[1]
  t1.expectKind(nnkBracketExpr)
  let t11 = t1[1]
  t11.expectKind(nnkInfix) # ..
  return t11[1]

macro rangeLow*(r: range):auto =
  #echo r.treerepr
  #echo r.getType.treerepr
  #echo r.getTypeImpl.treerepr
  #return newCall(ident("low"),r.getTypeImpl[1])
  let t = r.getTypeImpl
  t.expectKind(nnkBracketExpr)
  let t1 = t[1]
  t1.expectKind(nnkInfix) # ..
  return t1[1]

#template low*(r: typedesc[range]):untyped = r

#macro high*(r: typedesc[range]):auto =
#  echo r.treerepr
#  echo r.getType.treerepr
#  return r

macro rangeLen*(r: typedesc[range]):auto =
  #echo r.treerepr
  #echo r.getType.treerepr
  let t = r.getTypeImpl
  #echo t.treerepr
  t.expectKind(nnkBracketExpr)
  let t1 = t[1]
  t1.expectKind(nnkBracketExpr)
  let t11 = t1[1]
  t11.expectKind(nnkInfix) # ..
  return newCall(ident"-",newCall(ident"+",t11[2],newLit(1)),t11[1])

macro rangeLen*(r: range):auto =
  #echo r.treerepr
  #echo r.getType.treerepr
  let t = r.getTypeImpl
  #echo t.treerepr
  t.expectKind(nnkBracketExpr)
  let t1 = t[1]
  t1.expectKind(nnkInfix) # ..
  return newCall(ident"-",newCall(ident"+",t1[2],newLit(1)),t1[1])

macro echoVars*(vars:varargs[typed]):untyped =
  result = newStmtList()
  let e = ident "echo"
  for v in vars:
    result.add newCall(e, newLit($v & ": "), v)

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
