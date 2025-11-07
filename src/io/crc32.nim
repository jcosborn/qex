#import strutils

type Crc32* = uint32
const InitCrc32* = Crc32(not 0'u32)
const crc32PolyLow = Crc32(0xedb88320)
#const crc32Poly = uint64(0x100000000) + uint64(crc32PolyLow)

proc createCrcTable(): array[0..255, Crc32] =
  for i in 0..255:
    var rem = Crc32(i)
    for j in 0..7:
      if (rem and 1) > 0'u32: rem = (rem shr 1) xor crc32PolyLow
      else: rem = rem shr 1
    result[i] = rem

const crc32table = createCrcTable()

proc updateCrc32*(crc: Crc32, c: char): Crc32 =
  (crc shr 8) xor crc32table[(crc and 0xff) xor uint32(ord(c))]

proc updateCrc32*(crc: Crc32, buf: pointer, bytes: int): Crc32 =
  let cbuf = cast[ptr UncheckedArray[char]](buf)
  result = crc
  for i in 0..<bytes:
    result = updateCrc32(result, cbuf[i])

proc finishCrc32*(c: Crc32): Crc32 =
  not c

proc crc32Raw*(s: string, init=InitCrc32): Crc32 =
  result = init
  for c in s:
    result = updateCrc32(result, c)

proc crc32*(s: string): Crc32 =
  ## Compute the Crc32 on the string `s`
  result = crc32Raw(s, InitCrc32)
  result = finishCrc32(result)

proc polyMul(x0: uint32, y0: uint64): uint64 =
  var x = x0
  var y = y0
  #echo x.toHex, "  ", y.toHex
  while x!=0:
    if (x and 1)!=0:
      result = result xor y
    #echo x.toHex, "  ", result.toHex
    y = y shl 1
    x = x shr 1

#[
proc polyRem(x0: uint64, y0: uint64): uint32 =
  var x = x0
  var y = y0
  var b = 1.uint32
  var q = 0.uint32
  while (y and 0x8000000000000000'u64) == 0:
    y = y shl 1
    b = b shl 1
  while b != 0:
    #echo b, "  ", x, "  ", y
    if (x xor y) < x:
      x = x xor y
      q = q xor b
    y = y shr 1
    b = b shr 1
  echo q.toHex
  echo (x0 xor polyMul(q, y0)).toHex
  echo x.toHex
  result = x.uint32
]#

proc mulRem(r1,r2: Crc32): Crc32 =
  var t = polyMul(r1, r2) shl 1
  for i in 0..<4:
    result = updateCrc32(result, char(t and 255))
    t = t shr 8
  result = result xor Crc32(t)

#[
proc zeroPadCrc32X(crc: Crc32, n: int): Crc32 =
  var fac = Crc32(0x80000000)
  for i in 0..<n:
    fac = updateCrc32(fac, '\0')
  #echo "fac: ", fac.toHex
  result = mulRem(fac, crc)
]#

proc zeroPadCrc32*(crc: Crc32, n: int): Crc32 =
  var fac = Crc32(0x80000000)
  var s = Crc32(0x00800000)
  var nn = n
  while nn > 0:
    if (nn and 1) != 0:
      fac = mulRem(fac, s)
    s = mulRem(s, s)
    nn = nn shr 1
  #echo "fac: ", fac.toHex
  result = mulRem(fac, crc)

when isMainModule:
  echo "initCrc32 = ", $InitCrc32
  let s = "The quick brown fox jumps over the lazy dog"
  let foo = crc32(s)
  echo foo
  doAssert(foo == 0x414FA339)

  for i in 0..<100:
    doAssert(zeroPadCrc32(1,i) == zeroPadCrc32X(1,i))

  var n = s.len
  var n2 = n div 2
  var s1 = s[0..<n2]
  var s2 = s[n2..<n]
  let foo1x = crc32Raw(s1, InitCrc32)
  let foo1 = zeroPadCrc32(foo1x, n-n2)
  let foo2 = crc32Raw(s2, 0)
  let foo12 = finishCrc32(foo1 xor foo2)
  doAssert(foo12 == foo)

  var f = newSeq[Crc32](n)
  f[0] = InitCrc32
  for i in 0..<n:
    f[i] = crc32Raw(s[i..i], f[i])
    let m = n-1-i
    for j in 0..<m:
      f[i] = zeroPadCrc32(f[i], 1)
  var ff = Crc32(0)
  for i in 0..<n:
    ff = ff xor f[i]
  ff = finishCrc32(ff)
  doAssert(ff == foo)
