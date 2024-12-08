import base
import layout
export PerfInfo

type Fat7lCoefs* = object
  oneLink*: float
  threeStaple*: float
  fiveStaple*: float
  sevenStaple*: float
  lepage*: float

#type Fat7lState*[I,O] = object
#  coefs*: Fat7lCoefs
#  in*: I
#  out*: O

proc `$`*(c: Fat7lCoefs): string =
  result  = "oneLink:     " & $c.oneLink & "\n"
  result &= "threeStaple: " & $c.threeStaple & "\n"
  result &= "fiveStaple:  " & $c.fiveStaple & "\n"
  result &= "sevenStaple: " & $c.sevenStaple & "\n"
  result &= "lepage:      " & $c.lepage & "\n"

proc computeGenStaple(staple: auto, mu,nu: int, link: auto, coef: float,
                      gauge: auto, fl: auto, ts0: auto, ts1: auto,
                      tmat1: auto, tmat2: auto, ts2: auto) =
  ## Computes the staple:
  ##               mu
  ##            +-------+
  ##        nu  |       |
  ##            |       |
  ##            X       X
  ##  Where the mu link can be any su3_matrix. The result is saved in staple.
  ##  if staple==NULL then the result is not saved.
  ##  It also adds the computed staple to the fatlink[mu] with weight coef.
  mixin adj
  # Upper staple
  if link!=gauge[mu]:
    #QDP_M_eq_sM(ts0, link, QDP_neighbor[nu], QDP_forward, QDP_all);
    discard ts0 ^* link

  when type(staple) isnot bool: # Save the staple
    #QDP_M_eq_M_times_Ma(tmat1, ts0, ts1, QDP_all);
    #QDP_M_eq_M_times_M(staple, gauge[nu], tmat1, QDP_all);
    staple := gauge[nu] * ts0.field * ts1.field.adj
  else:  # No need to save the staple. Add it to the fatlinks
    #QDP_M_eq_M_times_Ma(tmat1, ts0, ts1, QDP_all);
    #QDP_M_eq_M_times_M(tmat2, gauge[nu], tmat1, QDP_all);
    #QDP_M_peq_r_times_M(fl[mu], &coef, tmat2, QDP_all);
    fl[mu] += coef * gauge[nu] * ts0.field * ts1.field.adj

  # lower staple
  #QDP_M_eq_Ma_times_M(tmat1, gauge[nu], link, QDP_all);
  #QDP_M_eq_M_times_M(tmat2, tmat1, ts1, QDP_all);
  #QDP_M_eq_sM(ts2, tmat2, QDP_neighbor[nu], QDP_backward, QDP_all);
  tmat2 := gauge[nu].adj * link * ts1.field
  discard ts2 ^* tmat2

  when type(staple) isnot bool: # Save the staple
    #QDP_M_peq_M(staple, ts2, QDP_all);
    #QDP_M_peq_r_times_M(fl[mu], &coef, staple, QDP_all);
    staple += ts2.field
    fl[mu] += coef * staple
  else:  # No need to save the staple. Add it to the fatlinks
    #QDP_M_peq_r_times_M(fl[mu], &coef, ts2, QDP_all);
    fl[mu] += coef * ts2.field

  #if(link!=gauge[mu]) QDP_discard_M(ts0);
  #QDP_discard_M(ts2);

proc makeImpLinks*(fl: auto, gf: auto, coef: auto,
                   ll: auto, gfLong: auto, naik: auto,
                   info: var PerfInfo) =
  tic("makeImpLinks")
  type lcm = type(gf[0])
  proc QDP_create_M(): lcm = result.new(gf[0].l)
  var
    staple: lcm
    tempmat1: lcm
    t1: lcm
    t2: lcm
    tsl: array[4,Shifter[lcm,type(gf[0][0])]]
    tsg: array[4,array[4,Shifter[lcm,type(gf[0][0])]]]
    ts1: array[4,Shifter[lcm,type(gf[0][0])]]
    ts2: array[4,Shifter[lcm,type(gf[0][0])]]
    nflop = 0.0
    #dtime = 0.0
    coef1 = coef.oneLink
    coef3 = coef.threeStaple
    coef5 = coef.fiveStaple
    coef7 = coef.sevenStaple
    coefL = coef.lepage
    have5 = (coef5!=0.0) or (coef7!=0.0) or (coefL!=0.0)
    have3 = (coef3!=0.0) or have5

  # to fix up the Lepage term, included by a trick below
  coef1 -= 6.0*coefL

  if have3 or naik!=0.0:
    nflop = 61632
    staple = QDP_create_M()
    tempmat1 = QDP_create_M()
    #if have3:
    if true:
      t1 = QDP_create_M()
      t2 = QDP_create_M()
      for dir in 0..<4:
        tsl[dir] = newShifter(gf[dir], dir, 1)
        ts1[dir] = newShifter(gf[dir], dir, 1)
        ts2[dir] = newShifter(gf[dir], dir, -1)
        for nu in 0..<4:
          if dir!=nu:
            tsg[dir][nu] = newShifter(gf[dir], nu, 1)
      threads:
        for dir in 0..<4:
          for nu in 0..<4:
            if dir!=nu:
              discard tsg[dir][nu] ^*! gf[dir]

  toc("main loop")
  threads:
    for dir in 0..<4:
      #QDP_M_eq_r_times_M(fl[dir], &coef1, gf[dir], QDP_all);
      fl[dir] := coef1 * gf[dir]
      if have3:
        for nu in 0..<4:
          if nu!=dir:
            compute_gen_staple(staple, dir, nu, gf[dir], coef3, gf, fl,
                               tsg[dir][nu], tsg[nu][dir], t1, t2, ts2[nu])
            if coefL!=0.0:
              compute_gen_staple(false, dir, nu, staple, coefL, gf, fl,
                                 tsl[nu], tsg[nu][dir], t1, t2, ts2[nu])
            if coef5!=0.0 or coef7!=0.0:
              for rho in 0..<4:
                if (rho!=dir) and (rho!=nu):
                  compute_gen_staple(tempmat1, dir, rho, staple, coef5, gf, fl,
                                     tsl[rho], tsg[rho][dir], t1, t2, ts2[rho])
                  if coef7!=0.0:
                    for sig in 0..<4:
                      if (sig!=dir) and (sig!=nu) and (sig!=rho):
                        compute_gen_staple(false, dir, sig, tempmat1,coef7,gf,fl,
                                           ts1[sig],tsg[sig][dir],t1,t2,ts2[sig])

    # long links
    if naik!=0.0:
      for dir in 0..<4:
        #QDP_M_eq_sM(staple, gfLong[dir], QDP_neighbor[dir], QDP_forward,QDP_all)
        #QDP_M_eq_M_times_M(tempmat1, gfLong[dir], staple, QDP_all)
        #QDP_M_eq_sM(staple, tempmat1, QDP_neighbor[dir], QDP_forward, QDP_all)
        #QDP_M_eq_M_times_M(ll[dir], gfLong[dir], staple, QDP_all)
        #QDP_M_eq_r_times_M(ll[dir], &naik, ll[dir], QDP_all)
        discard tsl[dir] ^* gfLong[dir]
        discard ts1[dir] ^* (gfLong[dir] * tsl[dir].field)
        ll[dir] := naik * (gfLong[dir] * ts1[dir].field)

  toc("end")
  inc info.count
  info.flops += nflop * gf[0].l.localGeom.prod
  info.secs += getElapsedTime()

proc makeImpLinks*(fl: auto, gf: auto, coef: auto, info: var PerfInfo) =
  makeImpLinks(fl, gf, coef, fl, gf, 0.0, info)

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
  var coef: Fat7lCoefs
  coef.oneLink = 1.0
  coef.threeStaple = 1.0
  coef.fiveStaple = 1.0
  coef.sevenStaple = 1.0
  coef.lepage = 1.0
  var naik = 1.0

  var fl = lo.newGauge()
  var ll = lo.newGauge()
  var g3 = g

  echo g.plaq
  makeImpLinks(fl, g, coef, info)
  makeImpLinks(fl, g, coef, ll, g3, naik, info)
  echo fl.plaq
  echo ll.plaq
  echo pow(1.0,4)/6.0
  echo pow(1.0+6.0,4)/6.0
  echo pow(1.0+6.0+6.0*4.0,4)/6.0
  echo pow(1.0+6.0+6.0*4.0+6.0*4.0*2.0,4)/6.0
  echo pow(1.0+6.0+6.0*4.0+6.0*4.0*2.0+6.0,4)/6.0

  qexFinalize()
