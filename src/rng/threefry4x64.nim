#[
Threefry4x64/20
John K. Salmon, Mark A. Moraes, Ron O. Dror, David E. Shaw
Parallel Random Numbers: As Easy as 1, 2, 3
Proceedings of SC11, 2011.
]#

import math
import comms/comms

type
  State = array[4,uint64]

  Threefry4x64* = object
    ## The 20-round Threefry4x64 counter generator as a QEX RNG.
    ## Each four-word block supplies eight uint32 values to `next`.
    c,k,o: State
    i: uint8

template mix2(x0,x1:var uint64; rx:SomeInteger;
              z0,z1:var uint64; rz:SomeInteger):auto =
  const
    x = 64 - rx
    z = 64 - rz
  x0 += x1
  z0 += z1
  x1 = (x1 shl rx) or (x1 shr x)
  z1 = (z1 shl rz) or (z1 shr z)
  x1 = x1 xor x0
  z1 = z1 xor z0

template mixk(x0,x1:var uint64; rx:SomeInteger;
              z0,z1:var uint64; rz:SomeInteger;
              k0,k1,l0,l1:uint64):auto =
  const
    x = 64 - rx
    z = 64 - rz
  x1 += k1
  z1 += l1
  x0 += x1+k0
  z0 += z1+l0
  x1 = (x1 shl rx) or (x1 shr x)
  z1 = (z1 shl rz) or (z1 shr z)
  x1 = x1 xor x0
  z1 = z1 xor z0

proc encrypt(r:var Threefry4x64) =
  var
    b0 = r.c[0]
    b1 = r.c[1]
    b2 = r.c[2]
    b3 = r.c[3]
  let
    k0 = r.k[0]
    k1 = r.k[1]
    k2 = r.k[2]
    k3 = r.k[3]
    k4 = 0x1BD11BDAA9FC1A22u64 xor k0 xor k1 xor k2 xor k3
  mixk(b0, b1, 14,   b2, b3, 16,   k0, k1, k2, k3)
  mix2(b0, b3, 52,   b2, b1, 57)
  mix2(b0, b1, 23,   b2, b3, 40)
  mix2(b0, b3,  5,   b2, b1, 37)
  mixk(b0, b1, 25,   b2, b3, 33,   k1, k2, k3, k4+1)
  mix2(b0, b3, 46,   b2, b1, 12)
  mix2(b0, b1, 58,   b2, b3, 22)
  mix2(b0, b3, 32,   b2, b1, 32)

  mixk(b0, b1, 14,   b2, b3, 16,   k2, k3, k4, k0+2)
  mix2(b0, b3, 52,   b2, b1, 57)
  mix2(b0, b1, 23,   b2, b3, 40)
  mix2(b0, b3,  5,   b2, b1, 37)
  mixk(b0, b1, 25,   b2, b3, 33,   k3, k4, k0, k1+3)

  mix2(b0, b3, 46,   b2, b1, 12)
  mix2(b0, b1, 58,   b2, b3, 22)
  mix2(b0, b3, 32,   b2, b1, 32)

  mixk(b0, b1, 14,   b2, b3, 16,   k4, k0, k1, k2+4)
  mix2(b0, b3, 52,   b2, b1, 57)
  mix2(b0, b1, 23,   b2, b3, 40)
  mix2(b0, b3,  5,   b2, b1, 37)
  r.o[0] = b0 + k0
  r.o[1] = b1 + k1
  r.o[2] = b2 + k2
  r.o[3] = b3 + k3 + 5

proc incCounter(r:var Threefry4x64) =
  r.c[0].inc
  if r.c[0] != 0: return
  r.c[1].inc
  if r.c[1] != 0: return
  r.c[2].inc
  if r.c[2] != 0: return
  r.c[3].inc

proc incCounter(r:var Threefry4x64, z:uint64) =
  if z > not r.c[0]:
    r.c[1].inc
    if r.c[1] == 0:
      r.c[2].inc
      if r.c[2] == 0:
        r.c[3].inc
  r.c[0] += z

proc incCounter(r:var Threefry4x64, b:int) =
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

template maxInt*(x: typedesc[Threefry4x64]): int = int uint32.high
template high*(x: typedesc[Threefry4x64]): uint = uint uint32.high
template high*(x: Threefry4x64): uint = uint uint32.high
template numInts*(x: typedesc[Threefry4x64]): int = int 4294967296

template isWrapper*(x: Threefry4x64): bool = false
template isWrapper*(x: typedesc[Threefry4x64]): bool = false
template has*(x: typedesc[Threefry4x64], y: typedesc): bool = y is Threefry4x64

proc `$`*(x:Threefry4x64):string =
  "Threefry4x64(" & $x.c & " " & $x.k & " " & $x.o & " " & $x.i & ")"

proc skip*(prn:var Threefry4x64, offset:uint64, base=0) =
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

proc seedX(prn:var Threefry4x64, seed,subsequence:uint64) =
  # Put the seed in the low key word and separate subsequences by 2^76 calls.
  for i in 0..<4:
    prn.c[i] = 0
    prn.k[i] = 0
  prn.k[0] = seed
  prn.i = 0
  prn.encrypt
  prn.skip(subsequence, subsequenceBase)

proc seedIndep*(prn:var Threefry4x64; sed,index:auto) =
  seedX(prn, sed.uint64, index.uint64)

proc seed*(prn:var Threefry4x64; sed,index:auto) =
  ## The seed `sed` is broadcasted from rank 0.
  ## For independent seeding, use `seedIndep`.
  var ss = sed
  defaultComm.broadcast(ss.addr, sizeof(ss))
  seedIndep(prn, ss, index)

# Threefry returns four uint64 words per counter.  Return each word as its
# low uint32 followed by its high uint32.
proc nextI(prn:var Threefry4x64):uint32 {.inline.} =
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

proc next64(prn:var Threefry4x64):uint64 {.inline.} =
  let
    lo = prn.nextI.uint64
    hi = prn.nextI.uint64
  result = lo or (hi shl 32)

proc integer*(prn:var Threefry4x64):int =
  ## Return random integer from 0 to maxInt
  result = int prn.nextI

proc next*(prn:var Threefry4x64):uint =
  ## Return random integer from 0 to maxInt
  result = uint prn.nextI

proc uniform*(prn:var Threefry4x64):float =
  ## Return random number uniform on (0,1]
  ## Use Random123's `u01<double>(uint64_t)` mapping.
  result = float(prn.next64) * norm + 0.5 * norm

proc gaussian*(prn:var Threefry4x64):float =
  ## Gaussian normal deviate
  ## Probability distribution exp( -x\*x/2 ), so < x^2 > = 1
  let
    v = prn.uniform
    p = prn.uniform * 2.0 * PI
    r = sqrt(-2.0 * ln(v))
  result = r * cos(p)

import maths/types
# Only needed for non-vectorized RNGs.
template gaussian*(x:var auto, r:MaskedObj[Threefry4x64]) =
  mixin gaussian
  gaussian(x, r[])
