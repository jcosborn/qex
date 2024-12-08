import base
import layout
import gauge
#import strUtils
import fat7l
import smearutil

export PerfInfo

const keepProj {.boolDefine.} = true
when keepProj:
  static: echo "hypsmear: keeping projected fields"
else:
  static: echo "hypsmear: NOT keeping projected fields"

type HypCoefs* = object
  alpha1*: float
  alpha2*: float
  alpha3*: float

proc `$`*(c: HypCoefs): string =
  result = "Hyp{\n"
  result &= "  alpha1: " & $c.alpha1 & "\n"
  result &= "  alpha2: " & $c.alpha2 & "\n"
  result &= "  alpha3: " & $c.alpha3 & "\n"
  result &= "}"

template proj(x: auto) =
  for e in x:
    x[e].projectU x[e]

template proj(r: auto, x: auto) =
  for e in r:
    r[e].projectU x[e]

template projDeriv(r: auto, x: auto, c: auto) =
  for i in r:
    r[i].projectUDeriv(x[i], c[i])

template projDeriv(r: auto, u, x: auto, c: auto) =
  for i in r:
    r[i].projectUDeriv(u[i], x[i], c[i])

# L[mu][nu] = P( (1-a1)*g[mu] + 0.5*a1 SS(g[nu],g[mu]) )
# L2[mu][nu] = P( (1-a2)*g[mu] + 0.25*a2 sum{a,b!=mu,nu} SS(L[a][b],L[mu][b]) )
# fl[mu] = P( (1-a3)*g[mu] + a3/6 sum{nu!=mu} SS(L2[nu][mu],L2[mu][nu]) )
#proc smear*(coef: HypCoefs, gf: auto, fl: auto, ht: HypTemps,
#            info: var PerfInfo) =
proc smearGetForce*[G](coef: HypCoefs, gf: G, fl: G,
            info: var PerfInfo):auto =
  ## Note that the resulting proc, smearedForce, holds a reference to the input gauge gf.
  ## The correctness of the algorithm depends on gf remaining the same.
  ## On the contrary, any changes to the smeared gauge fl would have no effects to the force calculation.
  tic()
  type lcm = type(gf[0])
  let lo = gf[0].l
  proc newlcm: lcm = result.new(gf[0].l)
  var
    l1x = newFieldArray2(lo,lcm,[4,4],mu!=nu)
    l2x = newOneOf(l1x)
    flx = newFieldArray(lo,lcm,4)
    tm1: lcm
    sm1: array[4,Shifter[lcm,type(gf[0][0])]]
    s1: array[4,array[4,Shifter[lcm,type(gf[0][0])]]]
    #nflop = 61632.0
    #dtime = 0.0
  when keepProj:
    var
      l1 = newOneOf(l1x)
      l2 = newOneOf(l1x)
  else:
    var
      lp1 = newlcm()
      lp2 = newlcm()

  tm1 = newlcm()
  for mu in 0..<4:
    sm1[mu] = newShifter(gf[mu], mu, -1)
    for nu in 0..<4:
      if nu!=mu:
        s1[mu][nu] = newShifter(gf[mu], nu, 1)
  threads:
    for mu in 0..<4:
      for nu in 0..<4:
        if nu!=mu:
          discard s1[mu][nu] ^* gf[mu]

  let
    alp1 = coef.alpha1 / 2.0
    alp2 = coef.alpha2 / 4.0
    alp3 = coef.alpha3 / 6.0
    ma1 = 1 - coef.alpha1
    ma2 = 1 - coef.alpha2
    ma3 = 1 - coef.alpha3

  toc("prep")
  threads:
    for mu in 0..<4:
      for nu in 0..<4:
        if nu!=mu:
          l1x[mu,nu] := ma1 * gf[mu]
          symStaple(l1x[mu,nu], alp1, gf[nu], gf[mu],
                    s1[nu][mu], s1[mu][nu], tm1, sm1[nu])
          when keepProj:
            l1[mu,nu].proj l1x[mu,nu]
    toc("1")

    for mu in 0..<4:
      for nu in 0..<4:
        if nu!=mu:
          l2x[mu,nu] := ma2 * gf[mu]
          for a in 0..<4:
            if a!=mu and a!=nu:
              let b = 1+2+3-mu-nu-a
              when keepProj:
                template lp1:untyped = l1[a,b]
                template lp2:untyped = l1[mu,b]
              else:
                lp1.proj l1x[a,b]
                lp2.proj l1x[mu,b]
              discard s1[nu][mu] ^* lp1
              discard s1[mu][a] ^* lp2
              symStaple(l2x[mu,nu], alp2, lp1, lp2,
                        s1[nu][mu], s1[mu][a], tm1, sm1[a])
          when keepProj:
            l2[mu,nu].proj l2x[mu,nu]
    toc("2")

    for mu in 0..<4:
      flx[mu] := ma3 * gf[mu]
      for nu in 0..<4:
        if nu!=mu:
          when keepProj:
            template lp1:untyped = l2[nu,mu]
            template lp2:untyped = l2[mu,nu]
          else:
            lp1.proj l2x[nu,mu]
            lp2.proj l2x[mu,nu]
          discard s1[nu][mu] ^* lp1
          discard s1[mu][nu] ^* lp2
          symStaple(flx[mu], alp3, lp1, lp2,
                    s1[nu][mu], s1[mu][nu], tm1, sm1[nu])
      fl[mu].proj flx[mu]
  toc("threads end")

  proc smearedForce(f,chain:G) =
    tic("smearedF")
    # fₓₚₜ ← chainₘₖₕ d/dUₓₚₜ^*[Vₘₖₕ(U)^*] + chainₘₕₖ^* d/dUₓₚₜ^*[Vₘₕₖ(U)]
    var
      fl1 = newFieldArray2(lo,lcm,[4,4],mu!=nu)
      fl2 = newOneOf(fl1)
      fc = newFieldArray(lo,lcm,4)
      fs: array[4,Shifter[lcm,type(gf[0][0])]]
      tm2: lcm
    tm2 = newlcm()
    for mu in 0..<4:
      fs[mu] = newShifter(fc[mu], mu, 1)
    toc("prep")

    threads:
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            fl1[mu,nu] := 0
            fl2[mu,nu] := 0

      # proj flx → fl, fc ← chain
      for mu in 0..<4:
        fc[mu].projDeriv(flx[mu], chain[mu])
      # link (gf, l2) → flx, (f, fl2) ← fc
      for mu in 0..<4:
        f[mu] := ma3 * fc[mu]
        fc[mu] *= alp3
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            when keepProj:
              template lp1:untyped = l2[nu,mu]
              template lp2:untyped = l2[mu,nu]
            else:
              lp1.proj l2x[nu,mu]
              lp2.proj l2x[mu,nu]
            discard s1[nu][mu] ^* lp1
            discard s1[mu][nu] ^* lp2
            discard fs[nu] ^* fc[mu]
            symStapleDeriv(fl2[nu,mu], fl2[mu,nu],
                           lp1, lp2, s1[nu][mu], s1[mu][nu],
                           fc[mu], fs[nu], tm1, tm2, sm1[nu], sm1[mu])
      toc("1")

      # proj l2x → l2, fl2 ← fl2
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            when keepProj:
              fl2[mu,nu].projDeriv(l2[mu,nu], l2x[mu,nu], fl2[mu,nu])
            else:
              fl2[mu,nu].projDeriv(l2x[mu,nu], fl2[mu,nu])
      # link (gf, l1) → l2x, (f, fl1) ← fl2
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            f[mu] += ma2 * fl2[mu,nu]
            fl2[mu,nu] *= alp2
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            for a in 0..<4:
              if a!=mu and a!=nu:
                let b = 1+2+3-mu-nu-a
                when keepProj:
                  template lp1:untyped = l1[a,b]
                  template lp2:untyped = l1[mu,b]
                else:
                  lp1.proj l1x[a,b]
                  lp2.proj l1x[mu,b]
                discard s1[nu][mu] ^* lp1
                discard s1[mu][a] ^* lp2
                discard fs[a] ^* fl2[mu,nu]
                symStapleDeriv(fl1[a,b], fl1[mu,b],
                               lp1, lp2, s1[nu][mu], s1[mu][a],
                               fl2[mu,nu], fs[a], tm1, tm2, sm1[a], sm1[mu])
      toc("2")

      # proj l1x → l1, fl1 ← fl1
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            when keepProj:
              fl1[mu,nu].projDeriv(l1[mu,nu], l1x[mu,nu], fl1[mu,nu])
            else:
              fl1[mu,nu].projDeriv(l1x[mu,nu], fl1[mu,nu])
            discard s1[mu][nu] ^* gf[mu]
      # link gf → l1, f ← fl1
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            f[mu] += ma1 * fl1[mu,nu]
            fl1[mu,nu] *= alp1
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            discard fs[nu] ^* fl1[mu,nu]
            symStapleDeriv(f[nu], f[mu],
                           gf[nu], gf[mu], s1[nu][mu], s1[mu][nu],
                           fl1[mu,nu], fs[nu], tm1, tm2, sm1[nu], sm1[mu])
    toc("end")

  toc("end")
  smearedForce

proc smearPriv[G](coef: HypCoefs, gf: G, fl: G, info: var PerfInfo) {.codegenDecl:
    "__attribute__((noinline)) $# $#$#".} =
  # Avoid inlining or other compiler optimizations
  # in order to guarantee the change of the stack pointer,
  # such that Nim's GC is able to collect the memory.
  {.emit: "asm (\"\");".}
  var f = coef.smearGetForce(gf, fl, info)
  f = nil
proc smear*[G](coef: HypCoefs, gf: G, fl: G, info: var PerfInfo) =
  ## Try our best to release memory here.
  ## Sometimes it still requires a GC after this function returns.
  coef.smearPriv(gf, fl, info)
  qexGC()

#proc smear*(c: HypCoefs, gf: auto, fl: auto, info: var PerfInfo) =
#  var t = newHypTemps(gf)
#  smear(c, gf, fl, t, info)

proc smear*(c: HypCoefs, g: auto, fl: auto) =
  var info: PerfInfo
  c.smear(g, fl, info)

# (d/dX') < C' F + F' C > /2
proc deriv*(coef: HypCoefs, gf: auto, fl: auto, info: var PerfInfo) =
  tic()
  type lcm = type(gf[0])
  proc newlcm: lcm = result.new(gf[0].l)
  var
    l1: array[4,array[4,lcm]]
    l2: array[4,array[4,lcm]]
    tm1: lcm
    sm1: array[4,Shifter[lcm,type(gf[0][0])]]
    s1: array[4,array[4,Shifter[lcm,type(gf[0][0])]]]
    nflop = 61632.0
    dtime = 0.0

  tm1 = newlcm()
  for mu in 0..<4:
    sm1[mu] = newShifter(gf[mu], mu, -1)
    for nu in 0..<4:
      if nu!=mu:
        l1[mu][nu] = newlcm()
        l2[mu][nu] = newlcm()
        s1[mu][nu] = newShifter(gf[mu], nu, 1)
        discard s1[mu][nu] ^* gf[mu]

  let alp1 = coef.alpha1 / 2.0
  for mu in 0..<4:
    #fl[mu] := 0
    for nu in 0..<4:
      if nu!=mu:
        l1[mu][nu] := (1-coef.alpha1) * gf[mu]
        symStaple(l1[mu][nu], alp1, gf[nu], gf[mu],
                  s1[nu][mu], s1[mu][nu], tm1, sm1[nu])
        proj l1[mu][nu]
        #fl[mu] += l1[mu][nu]

  let alp2 = coef.alpha2 / 4.0
  for mu in 0..<4:
    #fl[mu] := 0
    for nu in 0..<4:
      if nu!=mu:
        l2[mu][nu] := (1-coef.alpha2) * gf[mu]
        for a in 0..<4:
          if a!=mu and a!=nu:
            let b = 1+2+3-mu-nu-a
            discard s1[nu][mu] ^* l1[a][b]
            discard s1[mu][a] ^* l1[mu][b]
            symStaple(l2[mu][nu], alp2, l1[a][b], l1[mu][b],
                      s1[nu][mu], s1[mu][a], tm1, sm1[a])
        proj l2[mu][nu]
        #fl[mu] += l1[mu][nu]

  let alp3 = coef.alpha3 / 6.0
  for mu in 0..<4:
    fl[mu] := (1-coef.alpha3) * gf[mu]
    for nu in 0..<4:
      if nu!=mu:
        discard s1[nu][mu] ^* l2[nu][mu]
        discard s1[mu][nu] ^* l2[mu][nu]
        symStaple(fl[mu], alp3, l2[nu][mu], l2[mu][nu],
                  s1[nu][mu], s1[mu][nu], tm1, sm1[nu])
    proj fl[mu]

  toc()

#proc deriv*(c: HypCoefs, g: auto, fl: auto) =
#  var info: PerfInfo
#  c.deriv(g, fl, info)

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

  var info: PerfInfo
  var coef: HypCoefs
  coef.alpha1 = 0.4
  coef.alpha2 = 0.5
  coef.alpha3 = 0.5

  echo coef

  var fl = lo.newGauge()
  coef.smear(g, fl, info)
  resetTimers()

  template disp(g: typed) =
    let p = g.plaq
    let sp = 2.0*(p[0]+p[1]+p[2])
    let tp = 2.0*(p[3]+p[4]+p[5])
    echo p
    echo sp
    echo tp
    echo trace(g[0])

  disp g
  tic("main")
  coef.smear(g, fl, info)
  toc("smear")
  disp fl
  #echo pow(1.0,4)/6.0
  #echo pow(1.0+6.0,4)/6.0
  #echo pow(1.0+6.0+6.0*4.0,4)/6.0
  #echo pow(1.0+6.0+6.0*4.0+6.0*4.0*2.0,4)/6.0
  #echo pow(1.0+6.0+6.0*4.0+6.0*4.0*2.0+6.0,4)/6.0

  echoTimers()
  qexFinalize()
