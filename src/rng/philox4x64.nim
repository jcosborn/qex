#[
Philox4x64/10
John K. Salmon, Mark A. Moraes, Ron O. Dror, David E. Shaw
Parallel Random Numbers: As Easy as 1, 2, 3
Proceedings of SC11, 2011.
]#

import math
import comms/comms

type
  State = array[4,uint64]
  Key = array[2,uint64]

  Philox4x64* = object
    ## The 10-round Philox4x64 counter generator as a QEX RNG.
    ## Each four-word block supplies eight uint32 values to `next`.
    c,o: State
    k: Key
    i: uint8

const
  m0 = 0xD2E7470EE14C6C93u64
  m1 = 0xCA5A826395121157u64
  w0 = 0x9E3779B97F4A7C15u64
  w1 = 0xBB67AE8584CAA73Bu64

proc mulHi(a,b:uint64):uint64 {.inline.} =
  # Compile with `-d:philox4x64Portable` to force the portable path.
  when not defined(philox4x64Portable) and
      (defined(gcc) or defined(llvm_gcc) or defined(clang)):
    result = 0
    {.emit: """__uint128_t p = `a`; p *= `b`; `result` = p >> 64;""".}
  else:
    # High half of a uint64 product using Random123's 32-bit decomposition.
    const mask = 0xFFFFFFFFu64
    let
      lo = a * b
      ahi = a shr 32
      alo = a and mask
      bhi = b shr 32
      blo = b and mask
      ahbl = ahi * blo
      albh = alo * bhi
      s = (ahbl and mask) + (albh and mask)
    result = ahi*bhi + (ahbl shr 32) + (albh shr 32) + (s shr 32)
    if (lo shr 32) < (s and mask):
      result.inc

template mix(x0,x1,x2,x3:var uint64; k0,k1:uint64):auto =
  let
    hi0 = mulHi(m0, x0)
    hi1 = mulHi(m1, x2)
    lo0 = m0 * x0
    lo1 = m1 * x2
  x0 = hi1 xor x1 xor k0
  x1 = lo1
  x2 = hi0 xor x3 xor k1
  x3 = lo0

template mixk(x0,x1,x2,x3:var uint64; k0,k1:var uint64):auto =
  k0 += w0
  k1 += w1
  mix(x0, x1, x2, x3, k0, k1)

proc encrypt(r:var Philox4x64) =
  var
    b0 = r.c[0]
    b1 = r.c[1]
    b2 = r.c[2]
    b3 = r.c[3]
    k0 = r.k[0]
    k1 = r.k[1]
  mix(b0, b1, b2, b3, k0, k1)
  mixk(b0, b1, b2, b3, k0, k1)
  mixk(b0, b1, b2, b3, k0, k1)
  mixk(b0, b1, b2, b3, k0, k1)
  mixk(b0, b1, b2, b3, k0, k1)
  mixk(b0, b1, b2, b3, k0, k1)
  mixk(b0, b1, b2, b3, k0, k1)
  mixk(b0, b1, b2, b3, k0, k1)
  mixk(b0, b1, b2, b3, k0, k1)
  mixk(b0, b1, b2, b3, k0, k1)
  r.o = [b0, b1, b2, b3]

proc incCounter(r:var Philox4x64) =
  r.c[0].inc
  if r.c[0] != 0: return
  r.c[1].inc
  if r.c[1] != 0: return
  r.c[2].inc
  if r.c[2] != 0: return
  r.c[3].inc

proc incCounter(r:var Philox4x64, z:uint64) =
  if z > not r.c[0]:
    r.c[1].inc
    if r.c[1] == 0:
      r.c[2].inc
      if r.c[2] == 0:
        r.c[3].inc
  r.c[0] += z

proc incCounter(r:var Philox4x64, b:int) =
  var j = b shr 6
  let
    z = 1u64 shl (b and 63)
    x = r.c[j]
  r.c[j] += z
  if r.c[j] >= x: return
  inc j
  while j < 4:
    r.c[j].inc
    if r.c[j] != 0: return
    inc j

const
  norm = 1.0 / 18446744073709551616.0
  subsequenceBase = 76

template maxInt*(x: typedesc[Philox4x64]): int = int uint32.high
template high*(x: typedesc[Philox4x64]): uint = uint uint32.high
template high*(x: Philox4x64): uint = uint uint32.high
template numInts*(x: typedesc[Philox4x64]): int = int 4294967296

template isWrapper*(x: Philox4x64): bool = false
template isWrapper*(x: typedesc[Philox4x64]): bool = false
template has*(x: typedesc[Philox4x64], y: typedesc): bool = y is Philox4x64

proc `$`*(x:Philox4x64):string =
  "Philox4x64(" & $x.c & " " & $x.k & " " & $x.o & " " & $x.i & ")"

proc skip*(prn:var Philox4x64, offset:uint64, base=0) =
  ## Advance the state by `offset * 2^base` calls to `next`.
  ## `base` must be nonnegative; jumps wrap modulo the generator period.
  doAssert base >= 0, "base must be nonnegative"
  if offset == 0: return
  if base < 3:
    let
      q = (offset shl base) and 7
      t = prn.i.uint64 + q
      n = (offset shr (3-base)) + (t shr 3)
    prn.i = uint8(t and 7)
    if n != 0:
      prn.incCounter n
      prn.encrypt
    return

  var
    b = base - 3
    s = offset
    changed = false
  while s > 0 and b < 256:
    if (s and 1) != 0:
      prn.incCounter b
      changed = true
    s = s shr 1
    inc b
  if changed:
    prn.encrypt

proc seedX(prn:var Philox4x64, seed,subsequence:uint64) =
  # Put the seed in the low key word and separate subsequences by 2^76 calls.
  for i in 0..<4:
    prn.c[i] = 0
  for i in 0..<2:
    prn.k[i] = 0
  prn.k[0] = seed
  prn.i = 0
  prn.encrypt
  prn.skip(subsequence, subsequenceBase)

proc seedIndep*(prn:var Philox4x64; sed,index:auto) =
  seedX(prn, sed.uint64, index.uint64)

proc seed*(prn:var Philox4x64; sed,index:auto) =
  ## The seed `sed` is broadcasted from rank 0.
  ## For independent seeding, use `seedIndep`.
  var ss = sed
  defaultComm.broadcast(ss.addr, sizeof(ss))
  seedIndep(prn, ss, index)

# Philox returns four uint64 words per counter.  Return each word as its
# low uint32 followed by its high uint32.
proc nextI(prn:var Philox4x64):uint32 {.inline.} =
  let j = prn.i shr 1
  if (prn.i and 1) == 0:
    result = uint32(prn.o[j] and 0xFFFFFFFFu64)
  else:
    result = uint32(prn.o[j] shr 32)
  prn.i.inc
  if prn.i == 8:
    prn.i = 0
    prn.incCounter
    prn.encrypt

proc next64(prn:var Philox4x64):uint64 {.inline.} =
  let
    lo = prn.nextI.uint64
    hi = prn.nextI.uint64
  result = lo or (hi shl 32)

proc integer*(prn:var Philox4x64):int =
  ## Return random integer from 0 to maxInt
  result = int prn.nextI

proc next*(prn:var Philox4x64):uint =
  ## Return random integer from 0 to maxInt
  result = uint prn.nextI

proc uniform*(prn:var Philox4x64):float =
  ## Return random number uniform on (0,1]
  ## Use Random123's `u01<double>(uint64_t)` mapping.
  result = float(prn.next64) * norm + 0.5 * norm

proc gaussian*(prn:var Philox4x64):float =
  ## Gaussian normal deviate
  ## Probability distribution exp( -x\*x/2 ), so < x^2 > = 1
  let
    v = prn.uniform
    p = prn.uniform * 2.0 * PI
    r = sqrt(-2.0 * ln(v))
  result = r * cos(p)

import maths/types
# Only needed for non-vectorized RNGs.
template gaussian*(x:var auto, r:MaskedObj[Philox4x64]) =
  mixin gaussian
  gaussian(x, r[])
