import gauge, gauge/fat7l, maths, physics/stagD, strUtils

type
  AsqtadCoefs* = object
    fat7*: Fat7lCoefs
    naik*: float

proc setAsqtadFat7*(c: var Fat7lCoefs, f7lf,naik: float) =
  c.oneLink = (1.0+3.0*f7lf+naik)/8.0
  c.threeStaple = -1.0/16.0
  c.fiveStaple = 1.0/64.0
  c.sevenStaple = -1.0/384.0
  c.lepage = -f7lf/16.0

proc init*(c: var AsqtadCoefs) =
  var f7lf = 1.0
  var naik = 1.0
  c.fat7.setAsqtadFat7(f7lf, naik)
  c.naik = -naik/24.0

proc `$`*(c: AsqtadCoefs): string =
  let f1 = $c.fat7
  result  = "fat7:\n  " & f1.replace("\n","\n  ")
  result.removeSuffix("  ")
  result &= "naik: " & $c.naik & "\n"

proc smear*(c: AsqtadCoefs, g: any, fl,ll: any) =
  var info: PerfInfo
  makeImpLinks(fl, g, c.fat7, ll, g, c.naik, info)

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

  var coef: AsqtadCoefs
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
