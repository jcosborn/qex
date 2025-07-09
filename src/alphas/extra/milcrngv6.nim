## Brief: RNG from MILC version 6
## 
## Author: James C. Osborn
## Author: Xiao-Yong Jin
## 
## Modified: Curtis Taylor Peterson <curtistaylorpetersonwork@gmail.com>
## 
## Details:
## This is a copy of milcrng.nim for alphas. There are two major differences:
## 
## a.) Equivalent to running milcrng.nim with FUELcompat=1
## b.) Uses double-precision floats to better match MILC (ranstuff.c)
## 
## C language random number generator for parallel processors
## exclusive or of feedback shift register and integer congruence
## generator.  Use a different multiplier on each generator, and make sure
## that fsr is initialized differently on each generator.

import comms/qmp
import maths/types
import simd/simdWrap

import math

type MilcRngv6* = object
  r0, r1, r2, r3, r4, r5, r6: uint32
  icState, multiplier: uint32
  iset: int32
  gset: float

template isWrapper*(x: MilcRngv6): bool = false
template isWrapper*(x: typedesc[MilcRngv6]): bool = false
template has*(x: typedesc[MilcRngv6], y: typedesc): bool = y is MilcRngv6
template numberType*(x: MilcRngv6): typedesc = uint32
template numberType*(x: typedesc[MilcRngv6]): typedesc = uint32
template simdLength*(x: typedesc[MilcRngv6]): untyped = 1
template getNc*(x: MilcRngv6): untyped = 0
template getNs*(x: MilcRngv6): untyped = 0
template `:=`*(x: MilcRngv6, y: MilcRngv6) = 
  x = y
template `:=`*(x: MilcRngv6, y: Indexed) = 
  x := y[]
template `[]`*(x: MilcRngv6, y: Simd): untyped = x
template `[]=`*(x: MilcRngv6, y: Simd, z: typed) = 
  x := z

proc `$`*(x:MilcRngv6):string =
  result = "MilcRngv6 r:[ " & $x.r0
  result &= " " & $x.r1
  result &= " " & $x.r2
  result &= " " & $x.r3
  result &= " " & $x.r4
  result &= " " & $x.r5
  result &= " " & $x.r6
  result &= " ]  icState: " & $x.icState & " multiplier: " & $x.multiplier
  result &= " iset: " & $x.iset & " gset: " & $x.gset

const
  INDX1 = 69607'u32
  INDX2 = 8'u32
  ADDEND = 12345'u32
  MASK = 0x00FFFFFF'u32
  NUMINTS = 0x01000000'u32
  SCALE = 1.0'f32 / 0x01000000.float

template maxInt*(x: MilcRngv6): int = int(MASK)
template maxInt*(x: typedesc[MilcRngv6]): int = int(MASK)
template high*(x: MilcRngv6): uint = uint(MASK)
template high*(x: typedesc[MilcRngv6]): uint = uint(MASK)
template numInts*(x: MilcRngv6): int = int(NUMINTS)
template numInts*(x: typedesc[MilcRngv6]): int = int(NUMINTS)

# Try to strictly follow QLA_seed_random.c
proc seedX(prn: var MilcRngv6; seed0,index: int32) =
  const
    INDX1 = 69607'i32
    INDX2 = 8'i32
    ADDEND = 12345'i32
    MASK = 0x00FFFFFF'i32
  # Seed the generator
  # "index" selects which random number generator - which multiplier
  var seed = seed0
  template set(x: uint32) =
    seed = (INDX1 + INDX2 * index) * seed + ADDEND
    x = cast[uint32]((seed shr 8) and MASK)
  set(prn.r0)
  set(prn.r1)
  set(prn.r2)
  set(prn.r3)
  set(prn.r4)
  set(prn.r5)
  set(prn.r6)
  seed = (INDX1 + INDX2 * index) * seed + ADDEND
  prn.icState = seed.uint32
  prn.multiplier = uint32(100005'i32 + 8'i32 * index)
  prn.iset = 1
  prn.gset = 0

proc seedX(prn: var MilcRngv6; seed0,index: uint32) = 
  seedX(prn, seed0.int32, index.int32)

proc seedIndep*(prn: var MilcRngv6; sed,index: auto) =
  seedX(prn, sed.uint32, index.uint32)

proc seed*(prn: var MilcRngv6; sed,index: auto) =
  ## The seed `sed` is broadcasted from rank 0.
  ## For independent seeding, use `seedIndep`.
  var ss = sed
  QMP_broadcast(ss.addr, sizeof(ss).csize_t)
  seedIndep(prn, ss, index)

func nextI(prn: var MilcRngv6): uint32 {.inline.} =
  ## internal routine to return next value
  let t = (((prn.r5 shr 7) or (prn.r6 shl 17)) xor
      ((prn.r4 shr 1) or (prn.r5 shl 23))) and MASK
  prn.r6 = prn.r5
  prn.r5 = prn.r4
  prn.r4 = prn.r3
  prn.r3 = prn.r2
  prn.r2 = prn.r1
  prn.r1 = prn.r0
  prn.r0 = t
  let s = prn.ic_state * prn.multiplier + ADDEND
  prn.icState = s
  result = t xor ((s shr 8) and MASK)

func skip*(prn: var MilcRngv6, c = 1) =
  for i in 1..c:
    discard prn.nextI

proc integer*(prn: var MilcRngv6): int =
  ## Return random integer from 0 to maxInt
  result = int prn.nextI

proc next*(prn: var MilcRngv6): uint =
  ## Return random integer from 0 to maxInt
  result = uint prn.nextI

proc uniform*(prn: var MilcRngv6): float =
  ## Return random number uniform on [0,1)
  ## The choice of including endpoints may vary among different RNGs
  let i = prn.nextI
  result = SCALE * float(i)

proc agaussian*(prn: var MilcRngv6): float =
  ## agaussian normal deviate
  ## Probability distribution exp( -x*x/2 ), so < x^2 > = 1
  if prn.iset != 0:
    prn.iset = 0
    var v1,v2,rsq: float
    while true:
      v1 = float(prn.uniform)
      v2 = float(prn.uniform)
      v1 = 2.0 * v1 - 1.0
      v2 = 2.0 * v2 - 1.0
      rsq = v1 * v1 + v2 * v2
      if not ((rsq >= 1.0) or (rsq == 0.0)): break
    var fac = sqrt(-2.0 * ln(rsq) / rsq )
    prn.gset = v1 * fac
    result = v2 * fac
  else:
    prn.iset = 1
    result = prn.gset

# Only needed for non-vectorized RNGs
template agaussian*(x: var auto, r: MaskedObj[MilcRngv6]) =
  mixin agaussian
  agaussian(x, r[])

when isMainModule:
  import qex
  qexInit()

  var sed = intParam("seed", 987654321)
  var nu = intParam("nu", 10)
  var ng = intParam("ng", 10)
  installHelpParam()
  var s: MilcRngv6

  echo "uniform"
  s.seed(1, sed)
  if nu < 0:
    var n = 0
    while true:
      let x = s.uniform
      inc n
      if x==0.0 or x==1.0:
        echo n, " ", x
        #break
  else:
    for i in 1..nu:
      echo i, "\t", s.uniform

  echo "agaussian"
  s.seed(1, sed)
  if ng < 0:
    var n = 0
    while true:
      let x = s.agaussian
      inc n
      if x==0.0 or abs(x)>=10.0:
        echo n, " ", x
        #break
  else:
    for i in 1..ng:
      echo i, "\t", s.agaussian

  qexFinalize()