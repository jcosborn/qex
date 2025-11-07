#[
MRG32k3a
Pierre L'Ecuyer
Good Parameter Sets for Combined Multiple Recursive Random Number Generators
Operations Research, 47, 1 (1999), 159-164.
]#

import math
import comms/comms

type
  State = array[3,uint32]
  Trans = array[3,State]

  MRG32k3a* = object
    s1,s2: State

proc squaremod(a:Trans, m:uint64):auto =
  var x:Trans
  for i in 0..<3:
    var t:array[3,uint64]
    for k in 0..<3:
      let aik = a[i][k].uint64
      t[0] = t[0] + (aik*a[k][0].uint64 mod m)
      t[1] = t[1] + (aik*a[k][1].uint64 mod m)
      t[2] = t[2] + (aik*a[k][2].uint64 mod m)
    x[i][0] = uint32(t[0] mod m)
    x[i][1] = uint32(t[1] mod m)
    x[i][2] = uint32(t[2] mod m)
  return x

proc squaremod(a:Trans, p:static int, m:uint64):auto =
  var x:array[p,Trans]
  x[0] = a
  for i in 1..<p:
    x[i] = x[i-1].squaremod m
  return x

proc matvecmod(a:Trans, v:var State, m:uint64) =
  let
    v0 = v[0].uint64
    v1 = v[1].uint64
    v2 = v[2].uint64
  v[0] = uint32(((a[0][0].uint64*v0) mod m + (a[0][1].uint64*v1) mod m + (a[0][2].uint64*v2) mod m) mod m)
  v[1] = uint32(((a[1][0].uint64*v0) mod m + (a[1][1].uint64*v1) mod m + (a[1][2].uint64*v2) mod m) mod m)
  v[2] = uint32(((a[2][0].uint64*v0) mod m + (a[2][1].uint64*v1) mod m + (a[2][2].uint64*v2) mod m) mod m)

const
  norm = 2.328306549295728e-10
  m1 = 4294967087.0
  m2 = 4294944443.0
  a12 = 1403580.0
  a13n = 810728.0
  a21 = 527612.0
  a23n = 1370589.0
  defaultSEED = 12345U32
  subsequenceBase = 76
  # sequenceBase = 127

  a1:Trans = [
    [0u32,            1,           0],
    [0u32,            0,           1],
    [uint32(m1-a13n), uint32(a12), 0]]
  a2:Trans = [
    [0u32,            1, 0],
    [0u32,            0, 1],
    [uint32(m2-a23n), 0, uint32(a21)]]

  maxpower2 = 190
  a1sq = a1.squaremod(maxpower2, m1.uint64)
  a2sq = a2.squaremod(maxpower2, m2.uint64)

when a1sq[76]!=[[82758667u32, 1871391091u32, 4127413238u32], [3672831523u32, 69195019u32, 1871391091u32], [3672091415u32, 3528743235u32, 69195019u32]]:
  {.error:"a1sq[76] wrong!".}
when a2sq[76]!=[[1511326704u32, 3759209742u32, 1610795712u32], [4292754251u32, 1511326704u32, 3889917532u32], [3859662829u32, 4292754251u32, 3708466080u32]]:
  {.error:"a2sq[76] wrong!".}

#template maxInt*(x: MRG32k3a): int = int 4294967086
template maxInt*(x: typedesc[MRG32k3a]): int = int 4294967086
template high*(x: typedesc[MRG32k3a]): uint = uint 4294967086
template high*(x: MRG32k3a): uint = uint 4294967086
#template numInts*(x: MRG32k3a): int = 4294967087
template numInts*(x: typedesc[MRG32k3a]): int = int 4294967087

template isWrapper*(x: MRG32k3a): bool = false
template isWrapper*(x: typedesc[MRG32k3a]): bool = false
template has*(x: typedesc[MRG32k3a], y: typedesc): bool = y is MRG32k3a

proc `$`*(x:MRG32k3a):string =
  "MRG32k3a(" & $x.s1 & " " & $x.s2 & ")"

proc skip*(prn: var MRG32k3a, offset:uint64, base=0) =
  var
    i = 0
    s = offset
  while s>0:
    if (s and 1)!=0:
      matvecmod(a1sq[base+i], prn.s1, m1.uint64)
      matvecmod(a2sq[base+i], prn.s2, m2.uint64)
    s = s shr 1
    inc i

proc seedX(prn:var MRG32k3a, seed,subsequence:uint64) =
  if seed!=0:
    var d1 = defaultSEED.uint64 * uint64(uint32(seed) xor 0x55555555U)
    var d2 = defaultSEED.uint64 * uint64(uint32(seed shr 32) xor 0xAAAAAAAAU)
    prn.s1[0] = uint32(d1 mod m1.uint64)
    prn.s1[1] = uint32(d2 mod m1.uint64)
    prn.s1[2] = uint32(d1 mod m1.uint64)
    prn.s2[0] = uint32(d2 mod m2.uint64)
    prn.s2[1] = uint32(d1 mod m2.uint64)
    prn.s2[2] = uint32(d2 mod m2.uint64)
  else:
    prn.s1[0] = defaultSEED
    prn.s1[1] = defaultSEED
    prn.s1[2] = defaultSEED
    prn.s2[0] = defaultSEED
    prn.s2[1] = defaultSEED
    prn.s2[2] = defaultSEED
  prn.skip(subsequence, subsequenceBase)

proc seedIndep*(prn: var MRG32k3a; sed,index: auto) =
  seedX(prn, sed.uint64, index.uint64)
proc seed*(prn: var MRG32k3a; sed,index: auto) =
  ## The seed `sed` is broadcasted from rank 0.
  ## For independent seeding, use `seedIndep`.
  var ss = sed
  defaultComm.broadcast(ss.addr, sizeof(ss))
  seedIndep(prn, ss, index)

#[
proc next0(prn: var MRG32k3a): float {.inline.} =
  ## Return random integer uniform on [1,m1]
  var p1,p2:float
  p1 = a12 * prn.s1[1].float - a13n * prn.s1[0].float
  p1 = p1 mod m1
  if p1<0.0:
    p1 += m1
  prn.s1[0] = prn.s1[1]
  prn.s1[1] = prn.s1[2]
  prn.s1[2] = p1.uint32

  p2 = a21 * prn.s2[2].float - a23n * prn.s2[0].float
  p2 = p2 mod m2
  if p2<0.0:
    p2 += m2
  prn.s2[0] = prn.s2[1]
  prn.s2[1] = prn.s2[2]
  prn.s2[2] = p2.uint32

  if p1<=p2:
    result = p1 - p2 + m1
  else:
    result = p1 - p2
]#

proc nextI(prn: var MRG32k3a): int {.inline.} =
  ## Return random integer uniform on [1,m1]
  const
    a12i = int a12
    a13ni = int a13n
    a21i = int a21
    a23ni = int a23n
    m1i = int m1
    m2i = int m2
  var p1,p2: int
  p1 = a12i * prn.s1[1].int - a13ni * prn.s1[0].int
  p1 = p1 mod m1i
  if p1<0:
    p1 += m1i
  prn.s1[0] = prn.s1[1]
  prn.s1[1] = prn.s1[2]
  prn.s1[2] = p1.uint32

  p2 = a21i * prn.s2[2].int - a23ni * prn.s2[0].int
  p2 = p2 mod m2i
  if p2<0:
    p2 += m2i
  prn.s2[0] = prn.s2[1]
  prn.s2[1] = prn.s2[2]
  prn.s2[2] = p2.uint32

  if p1<=p2:
    result = p1 - p2 + m1i
  else:
    result = p1 - p2

proc integer*(prn: var MRG32k3a): int =
  ## Return random integer from 0 to maxInt
  result = int(prn.nextI) - 1

proc next*(prn: var MRG32k3a): uint =
  ## Return random integer from 0 to maxInt
  result = uint(prn.nextI - 1)

#[
proc uniform*(prn:var MRG32k3a): float =
  ## Return random number uniform on (0,1)
  var p1,p2:float
  p1 = a12 * prn.s1[1].float - a13n * prn.s1[0].float
  p1 = p1 mod m1
  if p1<0.0:
    p1 += m1
  prn.s1[0] = prn.s1[1]
  prn.s1[1] = prn.s1[2]
  prn.s1[2] = p1.uint32

  p2 = a21 * prn.s2[2].float - a23n * prn.s2[0].float
  p2 = p2 mod m2
  if p2<0.0:
    p2 += m2
  prn.s2[0] = prn.s2[1]
  prn.s2[1] = prn.s2[2]
  prn.s2[2] = p2.uint32

  if p1<=p2:
    result = (p1 - p2 + m1) * norm
  else:
    result = (p1 - p2) * norm
]#
proc uniform*(prn: var MRG32k3a): float =
  ## Return random number uniform on (0,1)
  result = norm * prn.nextI.float

proc gaussian*(prn: var MRG32k3a): float =
  ## Gaussian normal deviate
  ## Probability distribution exp( -x*x/2 ), so < x^2 > = 1
  var v,p,r: float
  v = prn.uniform
  p = prn.uniform * 2.0 * PI
  r = sqrt(-2.0 * ln(v))
  result = r * cos(p)

import maths/types
# Only needed for non-vectorized RNGs.
template gaussian*(x: var auto, r: MaskedObj[MRG32k3a]) =
  mixin gaussian
  gaussian(x, r[])

when isMainModule:
  echo a1sq[76]
  echo a2sq[76]
  var s: MRG32k3a
  commsInit()
  for k in 0..<2:
    echo "uniform sequence ",k
    s.seed(17u64^13, k)
    echo s
    for i in 1..10:
      echo i, "\t", s.uniform
  for k in 0..<2:
    echo "gaussian sequence ",k
    s.seed(1234, k)
    echo s
    for i in 1..10:
      echo i, "\t", s.gaussian
