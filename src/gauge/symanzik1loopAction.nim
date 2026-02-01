##  gauge action for Symanzik improved 1x1 + 1x2 + 1x1x1

import base, layout, maths, gaugeAction
export PerfInfo

type Seq2d[T] = object
    shape: array[2,int]
    data: seq[T]
proc newSeq2d[T](a,b: int): Seq2d[T] =
  let n = a*b
  result.shape = [a,b]
  result.data = newSeq[T](n)
proc map*(s: Seq2d; i,j: int): int = i*s.shape[1] + j
proc `[]`*[T](s: Seq2d[T]; i,j: SomeInteger): T =
  let k = s.map(i,j)
  s.data[k]
proc `[]=`*[T](s: var Seq2d[T]; i,j: SomeInteger, y: auto) =
  let k = s.map(i,j)
  s.data[k] = y

proc symanzik1loopAction*[T](coeffs: GaugeActionCoeffs;
                             links: openarray[T];
                             info: var PerfInfo): tuple[space,time:float] =
  tic("symanzik1loopAction")
  type
    GF = type links[0]
    GM = type links[0][0]
    Shft = Shifter[GF,GM]
  let
    nc = links[0][0].getNc
    EQMTM = (nc * nc * (8 * nc - 2))
    PEQMTM = (nc * nc * (8 * nc))
    fac = 1.0 / nc
    plaq = fac * coeffs.plaq
    rect = fac * coeffs.rect
    pgm = fac * coeffs.pgm
    adpl = fac * fac * coeffs.adjplaq
    lo = links[0].l
    nd = lo.nDim
    td = nd - 1
  var
    acts, actt = 0.0
    nflops = newThreadSingle(0)
    U = newSeq[GF](nd)
    Uf = newSeq2d[Shft](nd,nd)
  for mu in 0..<nd:
    U[mu] = links[mu]
    for nu in 0..<nd:
      if nu != mu:
        #Uf[mu,nu] = links[0].newOneOf
        #QDP_M_eq_sM(Uf[mu][nu], U[mu], neighbor[nu], QDP_forward, sub)
        Uf[mu,nu] = newShifter(U[mu], nu, 1)
  threads:
    for mu in 0..<nd:
      for nu in 0..<nd:
        if nu != mu:
          discard Uf[mu,nu] ^*! U[mu]

  var
    #UUf: array[nd, array[nd, ptr QDP_ColorMatrix]]
    #fstpl: array[nd, array[nd, ptr QDP_ColorMatrix]]
    #bstpl0: array[nd, array[nd, ptr QDP_ColorMatrix]]
    #bstpl: array[nd, array[nd, ptr QDP_ColorMatrix]]
    UUf = newSeq2d[GF](nd,nd)
    fstpl = newSeq2d[GF](nd,nd)
    bstpl0 = newSeq2d[GF](nd,nd)
    bstpl = newSeq2d[Shft](nd,nd)
    tc = lo.Complex
  for mu in 1..<nd:
    for nu in 0..<mu:
      if pgm != 0.0:
        UUf[mu,nu] = links[0].newOneOf
        UUf[nu,mu] = links[0].newOneOf
        fstpl[mu,nu] = links[0].newOneOf
        fstpl[nu,mu] = links[0].newOneOf
        bstpl0[mu,nu] = links[0].newOneOf
        bstpl0[nu,mu] = links[0].newOneOf
        bstpl[mu,nu] = newShifter(bstpl0[mu,nu], nu, -1)
        bstpl[nu,mu] = newShifter(bstpl0[nu,mu], mu, -1)
      else:
        if mu == 1:
          UUf[mu,nu] = links[0].newOneOf
          UUf[nu,mu] = links[0].newOneOf
          if rect != 0.0:
            fstpl[mu,nu] = links[0].newOneOf
            fstpl[nu,mu] = links[0].newOneOf
            bstpl0[mu,nu] = links[0].newOneOf
            bstpl0[nu,mu] = links[0].newOneOf
            bstpl[mu,nu] = newShifter(bstpl0[mu,nu], nu, -1)
            bstpl[nu,mu] = newShifter(bstpl0[nu,mu], mu, -1)
        else:
          UUf[mu,nu] = UUf[1,0]
          UUf[nu,mu] = UUf[0,1]
          if rect != 0.0:
            fstpl[mu,nu] = fstpl[1,0]
            fstpl[nu,mu] = fstpl[0,1]
            bstpl0[mu,nu] = bstpl0[1,0]
            bstpl0[nu,mu] = bstpl0[0,1]
            #bstpl[mu,nu] = bstpl[1,0]
            #bstpl[nu,mu] = bstpl[0,1]
            bstpl[mu,nu] = newShifter(bstpl0[mu,nu], nu, -1)
            bstpl[nu,mu] = newShifter(bstpl0[nu,mu], mu, -1)

  threads:
    var plaqs, plaqt, rects, rectt, pgms, pgmt, adpls, adplt = 0.0
    for mu in 1..<nd:
      for nu in 0..<mu:
        #QDP_M_eq_M_times_M(UUf[mu,nu], U[mu], Uf[nu,mu], sub)
        #QDP_M_eq_M_times_M(UUf[nu,mu], U[nu], Uf[mu,nu], sub)
        UUf[mu,nu] := U[mu] * Uf[nu,mu].field
        UUf[nu,mu] := U[nu] * Uf[mu,nu].field
        nflops += 2 * EQMTM
        if adpl != 0.0:
          #QDP_C_eq_M_dot_M(tc, UUf[mu,nu], UUf[nu,mu], sub)
          #tc := localdot(UUf[mu,nu], UUf[nu,mu])
          for e in lo:
            tc[e] := dot(UUf[mu,nu][e], UUf[nu,mu][e])
          #var z: QLA_Complex
          #var z: typeof(sum(tc))
          var zr = 0.0
          if plaq != 0.0:
            #QDP_c_eq_sum_C(addr(z), tc, sub)
            #z = sum(tc)
            zr = sum(tc).re
          #else:
          #  QLA_c_eq_r(z, 0)
          #var r: QLA_Real
          #QDP_r_eq_norm2_C(addr(r), tc, sub)
          let r = tc.norm2
          nflops += 8 * nc * nc - 2 + 2 + 4
          if mu == td:
            plaqt += zr
            adplt += r
          else:
            plaqs += zr
            adpls += r
          #QDP_destroy_C(tc)
          ## QOP_printf0("adpl[%i][%i] = %g\n", mu, nu, r);
          ## QOP_printf0("adpls: %g  adplt: %g\n", adpls, adplt);
        elif plaq != 0.0:
          #var t: QLA_Real
          #QDP_r_eq_re_M_dot_M(addr(t), UUf[mu,nu], UUf[nu,mu], sub)
          let t = redot(UUf[mu,nu], UUf[nu,mu])
          nflops += 4 * nc * nc
          if mu == td:
            plaqt += t
          else:
            plaqs += t
        if rect != 0.0 or pgm != 0.0:
          #QDP_M_eq_Ma_times_M(bstpl0[mu,nu], U[nu], UUf[mu,nu], sub)
          bstpl0[mu,nu] := U[nu].adj * UUf[mu,nu]
          #QDP_M_eq_sM(bstpl[mu,nu], bstpl0[mu,nu], neighbor[nu], QDP_backward, sub)
          discard bstpl[mu,nu] ^* bstpl0[mu,nu]
          #QDP_M_eq_Ma_times_M(bstpl0[nu,mu], U[mu], UUf[nu,mu], sub)
          bstpl0[nu,mu] := U[mu].adj * UUf[nu,mu]
          #QDP_M_eq_sM(bstpl[nu,mu], bstpl0[nu,mu], neighbor[mu], QDP_backward, sub)
          discard bstpl[nu,mu] ^* bstpl0[nu,mu]
          #QDP_M_eq_M_times_Ma(fstpl[mu,nu], UUf[nu,mu], Uf[nu,mu], sub)
          fstpl[mu,nu] := UUf[nu,mu] * Uf[nu,mu].field.adj
          #QDP_M_eq_M_times_Ma(fstpl[nu,mu], UUf[mu,nu], Uf[mu,nu], sub)
          fstpl[nu,mu] := UUf[mu,nu] * Uf[mu,nu].field.adj
          nflops += 4 * EQMTM
          if rect != 0.0:
            #var
            #  t = 0.0
            #  tr = 0.0
            #QDP_r_eq_re_M_dot_M(addr(t), bstpl[mu,nu], fstpl[mu,nu], sub)
            let t0 = redot(bstpl[mu,nu].field, fstpl[mu,nu])
            #inc(tr, t)
            #QDP_r_eq_re_M_dot_M(addr(t), bstpl[nu,mu], fstpl[nu,mu], sub)
            let t1 = redot(bstpl[nu,mu].field, fstpl[nu,mu])
            #inc(tr, t)
            nflops += 2 * 4 * nc * nc
            if mu == td:
              rectt += t0 + t1
            else:
              rects += t0 + t1

    template combineb(x, a, b, c: auto) =
      #var t: QLA_Real
      #QDP_M_eq_Ma_times_M(UUf[0,1], bstpl[a,c], bstpl[b,c], sub)
      UUf[0,1] := bstpl[a,c].field.adj * bstpl[b,c].field
      #QDP_M_eq_M_times_M(UUf[1,0], UUf[0,1], Uf[a,b], sub)
      UUf[1,0] := UUf[0,1] * Uf[a,b].field
      #QDP_r_eq_re_M_dot_M(addr(t), Uf[b,a], UUf[1,0], sub)
      let t = redot(Uf[b,a].field, UUf[1,0])
      x += t
      nflops += 2 * EQMTM + 4 * nc * nc

    template combineb2(x, a, b, c, d: auto) =
      #var t: QLA_Real
      #QDP_M_eq_Ma_times_M(UUf[0,1], bstpl[a,c], bstpl[b,c], sub)
      UUf[0,1] := bstpl[a,c].field.adj * bstpl[b,c].field
      #QDP_M_peq_Ma_times_M(UUf[0,1], bstpl[a,d], bstpl[b,d], sub)
      UUf[0,1] += bstpl[a,d].field.adj * bstpl[b,d].field
      #QDP_M_eq_M_times_M(UUf[1,0], UUf[0,1], Uf[a,b], sub)
      UUf[1,0] := UUf[0,1] * Uf[a,b].field
      #QDP_r_eq_re_M_dot_M(addr(t), Uf[b,a], UUf[1,0], sub)
      let t = redot(Uf[b,a].field, UUf[1,0])
      x += t
      nflops += 2 * EQMTM + PEQMTM + 4 * nc * nc

    template combinefb(x, a, b, c: auto) =
      #var t: QLA_Real
      #QDP_M_eq_Ma_times_M(UUf[0,1], fstpl[a,c], fstpl[b,c], sub)
      UUf[0,1] := fstpl[a,c].adj * fstpl[b,c]
      #QDP_M_peq_Ma_times_M(UUf[0,1], bstpl[a,c], bstpl[b,c], sub)
      UUf[0,1] += bstpl[a,c].field.adj * bstpl[b,c].field
      #QDP_M_eq_M_times_M(UUf[1,0], UUf[0,1], Uf[a,b], sub)
      UUf[1,0] := UUf[0,1] * Uf[a,b].field
      #QDP_r_eq_re_M_dot_M(addr(t), Uf[b,a], UUf[1,0], sub)
      let t = redot(Uf[b,a].field, UUf[1,0])
      x += t
      nflops += 2 * EQMTM + PEQMTM + 4 * nc * nc

    template combinefb2(x, a, b, c, d: auto) =
      #var t: QLA_Real
      #QDP_M_eq_Ma_times_M(UUf[0,1], fstpl[a,c], fstpl[b,c], sub)
      UUf[0,1] := fstpl[a,c].adj * fstpl[b,c]
      #QDP_M_peq_Ma_times_M(UUf[0,1], bstpl[a,c], bstpl[b,c], sub)
      UUf[0,1] += bstpl[a,c].field.adj * bstpl[b,c].field
      #QDP_M_peq_Ma_times_M(UUf[0,1], fstpl[a,d], fstpl[b,d], sub)
      UUf[0,1] += fstpl[a,d].adj * fstpl[b,d]
      #QDP_M_peq_Ma_times_M(UUf[0,1], bstpl[a,d], bstpl[b,d], sub)
      UUf[0,1] += bstpl[a,d].field.adj * bstpl[b,d].field
      #QDP_M_eq_M_times_M(UUf[1,0], UUf[0,1], Uf[a,b], sub)
      UUf[1,0] := UUf[0,1] * Uf[a,b].field
      #QDP_r_eq_re_M_dot_M(addr(t), Uf[b,a], UUf[1,0], sub)
      let t = redot(Uf[b,a].field, UUf[1,0])
      x += t
      nflops += 2 * EQMTM + 3 * PEQMTM + 4 * nc * nc

    if pgm != 0.0:
      #  FIXME: only works for nd=4
      if nd != 4:
        qexError("symanzik1loopAction with parallelogram only works for nDim == 4")
      combinefb(pgms, 0, 1, 2)
      combineb(pgms, 0, 2, 1)
      combineb(pgms, 1, 2, 0)
      #  rest
      # combinefb(pgmt,0,3,1);
      # combinefb(pgmt,0,3,2);
      combinefb2(pgmt, 0, 3, 1, 2)
      combineb(pgmt, 0, 1, 3)
      combineb(pgmt, 0, 2, 3)
      combinefb(pgmt, 1, 2, 3)
      # combineb(pgmt,1,3,0);
      # combineb(pgmt,1,3,2);
      combineb2(pgmt, 1, 3, 0, 2)
      # combineb(pgmt,2,3,0);
      # combineb(pgmt,2,3,1);
      combineb2(pgmt, 2, 3, 0, 1)

    threadSingle:
      acts = plaq * plaqs + rect * rects + pgm * pgms + adpl * adpls
      actt = plaq * plaqt + rect * rectt + pgm * pgmt + adpl * adplt

  let act0 = lo.physVol*(coeffs.plaq + 2*coeffs.rect + coeffs.adjplaq)
  let act1 = lo.physVol*(4*coeffs.pgm)
  result.space = 0.5*(nd-1)*(nd-2)*act0 + act1 - acts
  result.time = (nd-1)*act0 + (nd-1)*act1 - actt
  #info.final_sec = QOP_time() - dtime
  #info.final_flop = nflops * QDP_sites_on_node
  #info.status = QOP_SUCCESS

when isMainModule:
  import qex
  import physics/qcdTypes
  #import matrixFunctions
  qexInit()
  var defaultGaugeFile = "l88.scidac"
  #let defaultLat = @[2,2,2,2]
  let defaultLat = @[8,8,8,8]
  #let defaultLat = @[8,8,8]
  #let defaultLat = @[8,8]
  defaultSetup()
  #for mu in 0..<g.len: g[mu] := 1
  g.unit
  #g.random

  proc test(g:auto) =
    tic("test")
    echo "Test C_plaq = 1"
    var pl = plaq(g)
    echo "plaq:"
    echo pl
    echo pl.sum
    var gc = GaugeActionCoeffs(plaq:1.5)
    toc("plaq")
    var ga = gaugeAction.gaugeAction1(gc,g)
    toc("ga1")
    var ga2 = gaugeAction.gaugeAction2(gc,g)
    toc("ga2")
    var ga3 = gaugeAction.gaugeAction3(gc,g)
    toc("ga3")
    var gaA = gaugeAction.actionA(gc,g)
    toc("aA")
    echo "ga: ", ga, "\t", ga2, "\t", ga3, "\t", gaA
    var pi: PerfInfo
    let t = gc.symanzik1loopAction(g, pi)
    echo t
    toc("end")

  test(g)

  var nfail = 0
  proc check(gc,g:auto) =
    let nd = g[0].l.nDim
    let np = nd*(nd-1) div 2
    let npgm = 4*nd
    let a0 = g[0].l.physVol*(np*gc.plaq + 2*np*gc.rect + npgm*gc.pgm)
    #echo "a0: ", a0
    var a = a0 + gc.gaugeAction2(g)
    if gc.adjplaq != 0.0:
      var gc0 = GaugeActionCoeffs(adjplaq:gc.adjplaq)
      let ap = gc0.actionA(g)
      a += ap
    var pi: PerfInfo
    let b0 = gc.symanzik1loopAction(g, pi)
    let b = b0.space + b0.time
    let d = abs(a-b)
    echo "b0: ", b0
    echo d
    if d > 1e-6:
      inc nfail
      echo "a: ", a
      echo "b: ", b

  proc runcheck(g:auto) =
    template docheck(pl,rc,pg,ad:float) =
      check(GaugeActionCoeffs(plaq:pl,rect:rc,pgm:pg,adjplaq:ad),g)
    docheck(2.3,0,0,0)
    docheck(0,2.3,0,0)
    docheck(0,0,2.3,0)
    docheck(0,0,0,2.3)
    docheck(2.3,3.7,0,0)
    docheck(2.3,0,3.7,0)
    docheck(2.3,0,0,3.7)
    docheck(0,2.3,3.7,0)
    docheck(0,2.3,0,3.7)
    docheck(0,0,2.3,3.7)
    docheck(2.3,3.7,4.1,0)
    docheck(2.3,3.7,0,4.1)
    docheck(2.3,0,3.7,4.1)
    docheck(0,2.3,3.7,4.1)
    docheck(2.3,3.7,4.1,4.3)

  runcheck(g)
  g.random
  runcheck(g)

  if nfail > 0:
    echo "*** number of failures: ", nfail
    qexError("failed")
  else:
    echo "all tests passed"

  echoProf()
  qexFinalize()
