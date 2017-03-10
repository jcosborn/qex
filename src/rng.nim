import math

# RNG from MILC version 6
# C language random number generator for parallel processors
# exclusive or of feedback shift register and integer congruence
# generator.  Use a different multiplier on each generator, and make sure
# that fsr is initialized differently on each generator.

type RngMilc6* = object
  r0, r1, r2, r3, r4, r5, r6: uint32
  icState, multiplier: uint32

const
  INDX1 = 69607'u32
  INDX2 = 8'u32
  ADDEND = 12345'u32
  MASK = 0x00FFFFFF'u32
  SCALE = 1.0'f32 / 0x01000000.float32

proc seedX(prn: var RngMilc6; seed0,index: uint32) =
  ## Seed the generator
  ## "index" selects which random number generator - which multiplier
  var seed = seed0
  template set(x: untyped): untyped =
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
proc seed*(prn: var RngMilc6; sed,index: any) {.inline.} =
  seedX(prn, sed.uint32, index.uint32)

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
  when false:
    template iset = rs.addend
    template gset = rs.scale
    if iset != 0:
      iset = 0
      var
        v1: cdouble
        v2: cdouble
        rsq: cdouble
      while true:
        v1 = cast[cdouble](prn.uniform)
        v2 = cast[cdouble](prn.uniform)
        v1 = 2.0 * v1 - 1.0
        v2 = 2.0 * v2 - 1.0
        rsq = v1 * v1 + v2 * v2
        if not ((rsq >= 1.0) or (rsq == 0.0)): break
      var fac: cdouble = sqrt(-2.0 * ln(rsq) / rsq )
      gset = v1 * fac
      result = v2 * fac
    else:
      iset = 1
      result = gset
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
