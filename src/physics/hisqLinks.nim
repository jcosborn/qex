import gauge, gauge/fat7l, maths, strUtils

type
  HisqCoefs* = object
    fat7first*: Fat7lCoefs
    fat7second*: Fat7lCoefs
    naik*: float

proc setHisqFat7*(c: var Fat7lCoefs, f7lf,naik: float) =
  c.oneLink = (1.0+3.0*f7lf+naik)/8.0
  c.threeStaple = -1.0/16.0
  c.fiveStaple = 1.0/64.0
  c.sevenStaple = -1.0/384.0
  c.lepage = -f7lf/16.0

proc init*(c: var HisqCoefs) =
  var f7lf = 0.0
  var naik = 1.0
  c.fat7first.setHisqFat7(f7lf, 0.0)
  c.fat7second.setHisqFat7(2.0-f7lf, naik)
  c.naik = -naik/24.0

proc `$`*(c: HisqCoefs): string =
  let f1 = $c.fat7first
  let f2 = $c.fat7second
  result  = "fat7first:\n  " & f1.replace("\n","\n  ")
  result.removeSuffix("  ")
  result &= "fat7second:\n  " & f2.replace("\n","\n  ")
  result.removeSuffix("  ")
  result &= "naik: " & $c.naik & "\n"

proc smear*(c: HisqCoefs, g: auto, fl,ll: auto, t1,t2: auto) =
  mixin projectU
  var info: PerfInfo
  makeImpLinks(t1, g, c.fat7first, info)
  for mu in 0..<4:
    for i in t2[mu]:
      projectU(t2[mu][i], t1[mu][i])
  makeImpLinks(fl, t2, c.fat7second, ll, t2, c.naik, info)
proc smear*(c: HisqCoefs, g: auto, fl,ll: auto) =
  var t1 = fl.newOneOf
  var t2 = fl.newOneOf
  c.smear(g, fl, ll, t1, t2)

when isMainModule:
  import qex
  import physics/qcdTypes
  import gauge
  qexInit()
  #var defaultGaugeFile = "l88.scidac"
  let defaultLat = @[8,8,8,8]
  defaultSetup()
  #for mu in 0..<g.len: g[mu] := 1
  #g.random

  var coef: HisqCoefs
  g.setBC
  g.stagPhase
  coef.init()

  var fl = lo.newGauge()
  var ll = lo.newGauge()

  echo g.plaq
  coef.smear(g, fl, ll)
  echo fl.plaq
  echo ll.plaq
  echo pow(1.0,4)/6.0
  echo pow(1.0+6.0,4)/6.0
  echo pow(1.0+6.0+6.0*4.0,4)/6.0
  echo pow(1.0+6.0+6.0*4.0+6.0*4.0*2.0,4)/6.0
  echo pow(1.0+6.0+6.0*4.0+6.0*4.0*2.0+6.0,4)/6.0
