import fat7l, smearutil
import base/profile, layout

let STAPLE_FLOPS = 3*198+216+18
let SIDE_FORCE_FLOPS = 7*198+3*216+3*18

#[
proc staple(auto out, auto in0, auto link0, int mu, int nu) =
#define link     ftmp0[0]
#define linkmu   ftmp[0][mu]
#define in       ftmp0[1]
#define innu     ftmp[1][nu]
#define linkin   mtmp[0]
#define back     btmp0[0]
#define backnu   btmp[0][nu]
#define linkinnu mtmp[1]

  QDP_M_eq_M(link, link0, QDP_all);
  QDP_M_eq_sM(linkmu, link, QDP_neighbor[mu], QDP_forward, QDP_all);
  QDP_M_eq_M(in, in0, QDP_all);
  QDP_M_eq_sM(innu, in, QDP_neighbor[nu], QDP_forward, QDP_all);
  QDP_M_eq_Ma_times_M(linkin, link, in, QDP_all);
  QDP_M_eq_M_times_M(back, linkin, linkmu, QDP_all);
  QDP_M_eq_sM(backnu, back, QDP_neighbor[nu], QDP_backward, QDP_all);
  QDP_M_eq_M_times_M(linkinnu, link, innu, QDP_all);
  QDP_discard_M(innu);
  QDP_M_peq_M_times_Ma(out, linkinnu, linkmu, QDP_all);
  QDP_discard_M(linkmu);
  QDP_M_peq_M(out, backnu, QDP_all);
  QDP_discard_M(backnu);
#define STAPLE_FLOPS (3*198+216+18)

#undef link
#undef linkmu
#undef in
#undef innu
#undef linkin
#undef back
#undef backnu
#undef linkinnu
}
]#

#[
static void
side_force(QDP_ColorMatrix *force, QDP_ColorMatrix *bot0, QDP_ColorMatrix *side0,
           QDP_ColorMatrix *top0, int mu, int nu, QDP_ColorMatrix *stpl)
{
#define side     ftmp0[0]
#define sidemu   ftmp[0][mu]
#define top      ftmp0[1]
#define topnu    ftmp[1][nu]
#define bot      ftmp0[2]
#define botnu    ftmp[2][nu]
#define sidebot  mtmp[0]
#define sidetop  mtmp[1]
#define topnusidebot  btmp0[0]
#define fbmu          btmp[0][mu]
#define botnusidetop  btmp0[1]
#define fmbmu         btmp[1][mu]
#define sidebotsidemu btmp0[2]
#define stm           btmp[2][nu]
#define botnusidemu   mtmp[2]
#define botsidemu     mtmp[3]

  // force += bot * sidemu * topnu+
  // force -= bot-mu+ * side-mu * topnu-mu
  // -= top <-> bot
  // stpl += side * botnu * sidemu+
  // stpl += side-nu+ * bot-nu * sidemu-nu
  QDP_M_eq_M(side, side0, QDP_all);
  QDP_M_eq_sM(sidemu, side, QDP_neighbor[mu], QDP_forward, QDP_all);
  QDP_M_eq_M(top, top0, QDP_all);
  QDP_M_eq_sM(topnu, top, QDP_neighbor[nu], QDP_forward, QDP_all);
  QDP_M_eq_M(bot, bot0, QDP_all);
  QDP_M_eq_sM(botnu, bot, QDP_neighbor[nu], QDP_forward, QDP_all);
  QDP_M_eq_Ma_times_M(sidebot, side, bot, QDP_all);
  QDP_M_eq_Ma_times_M(sidetop, side, top, QDP_all);
  QDP_M_eq_Ma_times_M(topnusidebot, topnu, sidebot, QDP_all);
  QDP_M_eq_sM(fbmu, topnusidebot, QDP_neighbor[mu], QDP_backward, QDP_all);
  QDP_M_eq_Ma_times_M(botnusidetop, botnu, sidetop, QDP_all);
  QDP_M_eq_sM(fmbmu, botnusidetop, QDP_neighbor[mu], QDP_backward, QDP_all);
  QDP_M_eq_M_times_M(sidebotsidemu, sidebot, sidemu, QDP_all);
  QDP_M_eq_sM(stm, sidebotsidemu, QDP_neighbor[nu], QDP_backward, QDP_all);
  QDP_M_eq_M_times_Ma(botnusidemu, botnu, sidemu, QDP_all);
  QDP_discard_M(botnu);
  QDP_M_peq_M_times_M(stpl, side, botnusidemu, QDP_all);
 //QDP_M_meq_M_times_Ma(force, top, botnusidemu, QDP_all);
  QDP_M_peq_M_times_Ma(force, top, botnusidemu, QDP_all);
  QDP_M_eq_M_times_M(botsidemu, bot, sidemu, QDP_all);
  QDP_discard_M(sidemu);
  QDP_M_peq_M_times_Ma(force, botsidemu, topnu, QDP_all);
  QDP_discard_M(topnu);
  //QDP_M_meq_Ma(force, fbmu, QDP_all);
  QDP_M_peq_Ma(force, fbmu, QDP_all);
  QDP_discard_M(fbmu);
  QDP_M_peq_Ma(force, fmbmu, QDP_all);
  QDP_discard_M(fmbmu);
  QDP_M_peq_M(stpl, stm, QDP_all);
  QDP_discard_M(stm);
#define SIDE_FORCE_FLOPS (7*198+3*216+3*18)

#undef side
#undef sidemu
#undef top
#undef topnu
#undef bot
#undef botnu
#undef sidebot
#undef sidetop
#undef topnusidebot
#undef fbmu
#undef botnusidetop
#undef fmbmu
#undef sidebotsidemu
#undef stm
#undef botnusidemu
#undef botsidemu
}
]#

proc fat7lDeriv(deriv: auto, gauge: auto, mid: auto, coef: Fat7lCoefs,
                llderiv: auto, llgauge: auto, llmid: auto, naik: float,
                perf: var PerfInfo) =
  tic("fat7lDeriv")
  var nflops = ThreadSingle[int](0)
  let coefL = coef.lepage
  let coef1 = coef.oneLink - 6*coefL
  let coef3 = coef.threeStaple
  let coef5 = coef.fiveStaple
  let coef7 = coef.sevenStaple
  let have5 = (coef5!=0.0) or (coef7!=0.0) or (coefL!=0.0)
  let have3 = (coef3!=0.0) or have5
  type lcm = type(gauge[0])
  let lo = gauge[0].l
  proc newlcm: lcm = result.new(gauge[0].l)
  var
    stpl3: array[4,lcm]
    mid5: array[4,lcm]
    stpl5,mid3: lcm
    s1: array[4,array[4,Shifter[lcm,type(gauge[0][0])]]]
    s1a,s1b: array[4,Shifter[lcm,type(gauge[0][0])]]
    tm1,tm2: lcm
    sm1: array[4,Shifter[lcm,type(gauge[0][0])]]
  if have3:
    if have5:
      for mu in 0..3:
        stpl3[mu] = newlcm()
        mid5[mu] = newlcm()
        s1b[mu] = newShifter(gauge[0], mu, 1)
    stpl5 = newlcm()
    mid3 = newlcm()
    for mu in 0..3:
      for nu in 0..3:
        if nu!=mu:
          s1[mu][nu] = newShifter(gauge[mu], nu, 1)
      s1a[mu] = newShifter(gauge[0], mu, 1)
      sm1[mu] = newShifter(gauge[0], mu, -1)
    tm1 = newlcm()
    tm2 = newlcm()
    threads:
      for mu in 0..<4:
        for nu in 0..<4:
          if nu!=mu:
            discard s1[mu][nu] ^*! gauge[mu]
  threads:
    for sig in 0..3:
      if have5:
        for mu in 0..3:
          if mu==sig: continue
          #QDP_M_eq_zero(stpl3[mu], QDP_all);
          stpl3[mu] := 0
          #staple(stpl3[mu], gauge[sig], gauge[mu], sig, mu);
          symStaple(stpl3[mu], 1.0, gauge[mu], gauge[sig],
                    s1[mu][sig], s1[sig][mu], tm1, sm1[mu])
          #QDP_M_eq_zero(mid5[mu], QDP_all);
          mid5[mu] := 0
          nflops += STAPLE_FLOPS
        for rho in 0..3:
          if rho==sig: continue
          #QDP_M_eq_zero(stpl5, QDP_all);
          stpl5 := 0
          if coef7!=0.0:
            for mu in 0..3:
              if  mu==sig or mu==rho: continue
              for nu in 0..3:
                if nu==sig or nu==rho or nu==mu: continue
                #staple(stpl5, stpl3[mu], gauge[nu], sig, nu);
                discard s1a[nu] ^* stpl3[mu]
                symStaple(stpl5, 1.0, gauge[nu], stpl3[mu],
                          s1[nu][sig], s1a[nu], tm1, sm1[nu])
                nflops += STAPLE_FLOPS
            #QDP_M_eq_r_times_M(stpl5, &coef7, stpl5, QDP_all);
            stpl5 *= coef7
            nflops += 18
          #QDP_M_eq_zero(mid3, QDP_all);
          mid3 := 0
          if coefL!=0.0:
            #QDP_M_peq_r_times_M(stpl5, &coefL, stpl3[rho], QDP_all);
            stpl5 += coefL * stpl3[rho]
            nflops += 36
          if coefL!=0.0 or coef7!=0.0:
            #side_force(deriv[rho], mid[sig], gauge[rho], stpl5, sig, rho, mid3);
            discard s1a[rho] ^* stpl5
            discard s1b[rho] ^* mid[sig]
            symStapleDeriv(deriv[rho], mid3,
                           gauge[rho], stpl5, s1[rho][sig], s1a[rho],
                           mid[sig], s1b[rho], tm1, tm2, sm1[rho], sm1[sig])
            nflops += SIDE_FORCE_FLOPS
          if coefL!=0.0:
            #QDP_M_peq_r_times_M(mid5[rho], &coefL, mid3, QDP_all);
            mid5[rho] += coefL * mid3
            nflops += 36
          #//QDP_M_eqm_r_times_M(mid3, &coef7, mid3, QDP_all);
          #QDP_M_eq_r_times_M(mid3, &coef7, mid3, QDP_all);
          mid3 *= coef7
          #QDP_M_peq_r_times_M(mid3, &coef5, mid[sig], QDP_all);
          mid3 += coef5 * mid[sig]
          nflops += 18+36
          for mu in 0..3:
            if mu==sig or mu==rho: continue
            for nu in 0..3:
              if nu==sig or nu==rho or nu==mu: continue
              #side_force(deriv[mu], mid3, gauge[mu], stpl3[nu], sig, mu, mid5[nu]);
              discard s1a[mu] ^* stpl3[nu]
              discard s1b[mu] ^* mid3
              symStapleDeriv(deriv[mu], mid5[nu],
                             gauge[mu], stpl3[nu], s1[mu][sig], s1a[mu],
                             mid3, s1b[mu], tm1, tm2, sm1[mu], sm1[sig])
              nflops += SIDE_FORCE_FLOPS

      if have3:
        #QDP_M_eq_zero(mid3, QDP_all);
        mid3 := 0
        for mu in 0..3:
          if mu==sig: continue
          if have5:
            #//QDP_M_eq_r_times_M_minus_M(stpl5, &coef3, mid[sig], mid5[mu], QDP_all);
            #QDP_M_eq_r_times_M_plus_M(stpl5, &coef3, mid[sig], mid5[mu], QDP_all);
            stpl5 := coef3*mid[sig] + mid5[mu]
            nflops += 36
          else:
            #QDP_M_eq_r_times_M(stpl5, &coef3, mid[sig], QDP_all);
            stpl5 := coef3*mid[sig]
            nflops += 18
          #side_force(deriv[mu], stpl5, gauge[mu], gauge[sig], sig, mu, mid3);
          discard s1a[mu] ^* stpl5
          symStapleDeriv(deriv[mu], mid3,
                         gauge[mu], gauge[sig], s1[mu][sig], s1[sig][mu],
                         stpl5, s1a[mu], tm1, tm2, sm1[mu], sm1[sig])
          nflops += SIDE_FORCE_FLOPS
        #//QDP_M_meq_M(deriv[sig], mid3, QDP_all);
        #QDP_M_peq_M(deriv[sig], mid3, QDP_all);
        deriv[sig] += mid3
        nflops += 18
      if coef1!=0.0:
        #QDP_M_peq_r_times_M(deriv[sig], &coef1, mid[sig], QDP_all);
        deriv[sig] += coef1 * mid[sig]
        nflops += 36

    if coefL!=0.0:
      # fix up Lepage term
      let fixL = -coefL
      for mu in 0..3:
        #QDP_M_eq_zero(ftmp0[0], QDP_all);
        tm1 := 0
        for nu in 0..3:
          if nu==mu: continue
          #QDP_M_eq_Ma_times_M(btmp0[0], mid[nu], gauge[nu], QDP_all);
          tm2 := mid[nu].adj * gauge[nu]
          #QDP_M_eq_sM(btmp[0][nu], btmp0[0], QDP_neighbor[nu], QDP_backward, QDP_all);
          discard sm1[nu] ^* tm2
          #QDP_M_eq_M_times_Ma(stpl5, mid[nu], gauge[nu], QDP_all);
          stpl5 := mid[nu] * gauge[nu].adj
          #QDP_M_meq_M(stpl5, btmp[0][nu], QDP_all);
          stpl5 -= sm1[nu].field
          #QDP_discard_M(btmp[0][nu]);
          #QDP_M_peq_M(ftmp0[0], stpl5, QDP_all);
          #QDP_M_peq_Ma(ftmp0[0], stpl5, QDP_all);
          tm1 += stpl5 + stpl5.adj
        #QDP_M_eq_sM(ftmp[0][mu], ftmp0[0], QDP_neighbor[mu], QDP_forward, QDP_all);
        discard s1a[mu] ^* tm1
        #QDP_M_eq_M_times_M(stpl5, ftmp0[0], gauge[mu], QDP_all);
        stpl5 := tm1 * gauge[mu]
        #QDP_M_meq_M_times_M(stpl5, gauge[mu], ftmp[0][mu], QDP_all);
        stpl5 -= gauge[mu] * s1a[mu].field
        #QDP_discard_M(ftmp[0][mu]);
        #QDP_M_peq_r_times_M(deriv[mu], &fixL, stpl5, QDP_all);
        deriv[mu] += fixL * stpl5
      nflops += 4*(3*(2*198+3*18)+198+216+36)

    if naik!=0.0:
      for mu in 0..3:
        # force += mid * Umumu' * Umu'
        # force -= U-mu' * mid-mu * Umu'
        # force += U-mu' * U-mu-mu' * mid-mu-mu
        #QDP_M_eq_M(Uf, U, QDP_all);
        #QDP_M_eq_sM(Umu, Uf, QDP_neighbor[mu], QDP_forward, QDP_all);
        discard s1a[mu] ^* llgauge[mu]
        #QDP_M_eq_Ma_times_M(Umid, Uf, mid, QDP_all);
        tm1 := llgauge[mu].adj * llmid[mu]
        #QDP_M_eq_sM(Umidbmu, Umid, QDP_neighbor[mu], QDP_backward, QDP_all);
        discard sm1[mu] ^* tm1
        #QDP_M_eq_Ma_times_Ma(UmuU, Umu, Uf, QDP_all);
        tm2 := s1a[mu].field.adj * llgauge[mu].adj
        #QDP_M_eq_sM(UmuUs, UmuU, QDP_neighbor[mu], QDP_forward, QDP_all);
        discard s1b[mu] ^* tm2
        #QDP_M_eq_Ma_times_M(f3b, Uf, Umidbmu, QDP_all);
        llderiv[mu] += naik*(sm1[mu].field*s1a[mu].field.adj)
        tm1 := llgauge[mu].adj * sm1[mu].field
        #QDP_M_eq_sM(f3, f3b, QDP_neighbor[mu], QDP_backward, QDP_all);
        discard sm1[mu] ^* tm1
        #QDP_M_eq_M_times_M(f, mid, UmuUs, QDP_all);
        #QDP_discard_M(UmuUs);
        #QDP_M_peq_M_times_Ma(f, Umidbmu, Umu, QDP_all);
        #QDP_discard_M(Umidbmu);
        #QDP_discard_M(Umu);
        #QDP_M_peq_M(f, f3, QDP_all);
        #QDP_discard_M(f3);
        #QDP_M_peq_r_times_M(deriv[mu], &coefN, f, QDP_all);
        llderiv[mu] += naik*(llmid[mu]*s1b[mu].field+sm1[mu].field)

  #info->final_sec = QOP_time() - dtime;
  #info->final_flop = nflops*QDP_sites_on_node;
  #info->status = QOP_SUCCESS;
  toc("end")
  inc perf.count
  perf.flops += nflops * gauge[0].l.localGeom.prod
  perf.secs += getElapsedTime()

proc fat7lDeriv(deriv: auto, gauge: auto, mid: auto, coef: Fat7lCoefs,
                perf: var PerfInfo) =
  fat7lDeriv(deriv, gauge, mid, coef, deriv, gauge, mid, 0.0, perf)


when isMainModule:
  import qex, physics/qcdTypes
  import strformat
  qexInit()
  #var defaultGaugeFile = "l88.scidac"
  let defaultLat = @[8,8,8,8]
  defaultSetup()
  for mu in 0..<g.len: g[mu] := 1
  #g.random

  var info: PerfInfo
  var coef: Fat7lCoefs
  coef.oneLink = 1.0
  coef.threeStaple = 1.0
  coef.fiveStaple = 1.0
  coef.sevenStaple = 1.0
  coef.lepage = 1.0
  var naik = 1.0

  var fl = lo.newGauge()
  var fl2 = lo.newGauge()
  var ll = lo.newGauge()
  var ll2 = lo.newGauge()
  var dfl = lo.newGauge()
  var g2 = lo.newGauge()
  var dg = lo.newGauge()
  var fd = lo.newGauge()
  var ld = lo.newGauge()
  for mu in 0..<dg.len:
    dg[mu] := 0.00001 * g[mu]
    g2[mu] := g[mu] + dg[mu]
    fd[mu] := 0

  echo g.plaq
  echo g2.plaq
  makeImpLinks(fl, g, coef, info)
  info.clear
  resetTimers()
  makeImpLinks(fl, g, coef, info)
  echo info
  makeImpLinks(fl2, g2, coef, info)
  echo info
  echo fl.plaq
  echo fl2.plaq

  fat7lderiv(fd, g, dg, coef, info)
  echo info
  for mu in 0..3:
    dfl[mu] := fl2[mu] - fl[mu]
    #echo dfl[mu].norm2
    #echo fd[mu].norm2
    dfl[mu] -= fd[mu]
    echo dfl[mu].norm2

  for mu in 0..<dg.len:
    fd[mu] := 0
    ld[mu] := 0
  makeImpLinks(fl, g, coef, ll, g, naik, info)
  makeImpLinks(fl2, g2, coef, ll2, g2, naik, info)
  fat7lderiv(fd, g, dg, coef, ld, g, dg, naik, info)
  echo info
  for mu in 0..3:
    dfl[mu] := fl2[mu] - fl[mu]
    #echo dfl[mu].norm2
    #echo fd[mu].norm2
    dfl[mu] -= fd[mu]
    echo dfl[mu].norm2
    dfl[mu] := ll2[mu] - ll[mu]
    dfl[mu] -= ld[mu]
    echo dfl[mu].norm2

  echoProf()
  qexFinalize()
