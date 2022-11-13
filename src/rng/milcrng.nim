import math
import comms/qmp
import maths/types
import simd/simdWrap

# RNG from MILC version 6
# C language random number generator for parallel processors
# exclusive or of feedback shift register and integer congruence
# generator.  Use a different multiplier on each generator, and make sure
# that fsr is initialized differently on each generator.

type RngMilc6* = object
  r0, r1, r2, r3, r4, r5, r6: uint32
  icState, multiplier: uint32
  when defined(FUELCompat):
    # for Gaussian
    iset: int32
    gset: float

template isWrapper*(x: RngMilc6): untyped = false
template numberType*(x: RngMilc6): typedesc = uint32
template numberType*(x: typedesc[RngMilc6]): typedesc = uint32
template simdLength*(x: typedesc[RngMilc6]): untyped = 1
template getNc*(x: RngMilc6): untyped = 0
template getNs*(x: RngMilc6): untyped = 0
template `:=`*(x: RngMilc6, y: RngMilc6) =
  x = y
template `:=`*(x: RngMilc6, y: Indexed) =
  x := y[]
template `[]`*(x: RngMilc6, y: Simd): untyped = x
template `[]=`*(x: RngMilc6, y: Simd, z: typed) =
  x := z

#template `[]`*(x: RngMilc6): untyped = x

proc `$`*(x:RngMilc6):string =
  result = "RngMilc6 r:[ " & $x.r0
  result &= " " & $x.r1
  result &= " " & $x.r2
  result &= " " & $x.r3
  result &= " " & $x.r4
  result &= " " & $x.r5
  result &= " " & $x.r6
  result &= " ]  icState: " & $x.icState & " multiplier: " & $x.multiplier
  when defined(FUELCompat):
    result &= " iset: " & $x.iset & " gset: " & $x.gset

const
  INDX1 = 69607'u32
  INDX2 = 8'u32
  ADDEND = 12345'u32
  MASK = 0x00FFFFFF'u32
  SCALE = 1.0'f32 / 0x01000000.float32

when defined(FUELCompat):
  # Try to strictly follow QLA_seed_random.c
  proc seedX(prn: var RngMilc6; seed0,index: int32) =
    const
      INDX1 = 69607'i32
      INDX2 = 8'i32
      ADDEND = 12345'i32
      MASK = 0x00FFFFFF'i32
    ## Seed the generator
    ## "index" selects which random number generator - which multiplier
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
  proc seedX(prn: var RngMilc6; seed0,index: uint32) = seedX(prn, seed0.int32, index.int32)
else:
  proc seedX(prn: var RngMilc6; seed0,index: uint32) =
    ## Seed the generator
    ## "index" selects which random number generator - which multiplier
    var seed = seed0
    template set(x: uint32) =
      seed = (INDX1 + INDX2 * index) * seed + ADDEND
      x = (seed shr 8) and MASK
    set(prn.r0)
    set(prn.r1)
    set(prn.r2)
    set(prn.r3)
    set(prn.r4)
    set(prn.r5)
    set(prn.r6)
    seed = (INDX1 + INDX2 * index) * seed + ADDEND
    prn.icState = seed
    prn.multiplier = 100005'u32 + 8'u32 * index
    #prn.addend = 12345
    #prn.scale = 1.0 / float32(0x01000000)
proc seedIndep*(prn: var RngMilc6; sed,index: auto) {.inline.} =
  seedX(prn, sed.uint32, index.uint32)
proc seed*(prn: var RngMilc6; sed,index: auto) {.inline.} =
  ## The seed `sed` is broadcasted from rank 0.
  ## For independent seeding, use `seedIndep`.
  var ss = sed
  QMP_broadcast(ss.addr, sizeof(ss).csize_t)
  seedIndep(prn, ss, index)

proc uniform*(prn: var RngMilc6): float32 =
  ## Return random number uniform on [0,1]
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
  result = SCALE * (t xor ((s shr 8) and MASK)).float32

#var QLA_use_milc_gaussian* = false

proc gaussian*(prn: var RngMilc6): float32 =
  ## Gaussian normal deviate
  ## Probability distribution exp( -x*x/2 ), so < x^2 > = 1
  #if QLA_use_milc_gaussian:
  when defined(FUELCompat):
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
  else:
    const
      TINY = 9.999999999999999e-308
    var
      v: cdouble
      p: cdouble
      r: cdouble
    v = prn.uniform
    p = prn.uniform * 2.0 * PI
    r = sqrt(-2.0 * ln(v + TINY))
    result = r * cos(p)

# Only needed for non-vectorized RNGs.
template gaussian*(x: var auto, r: MaskedObj[RngMilc6]) =
  mixin gaussian
  gaussian(x, r[])

when isMainModule:
  var s: RngMilc6
  s.seed(1, 987654321)
  echo "uniform"
  for i in 1..10:
    echo i, "\t", s.uniform
  s.seed(1, 0)
  echo "gaussian"
  for i in 1..10:
    echo i, "\t", s.gaussian
