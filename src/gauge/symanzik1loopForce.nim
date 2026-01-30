##  gauge force for Symanzik improved 1x1 + 1x2 + 1x1x1

import base, layout, maths, gaugeAction, symanzik1loopAction
export PerfInfo

converter toField[T,U](x: Shifter[T,U]): T = x.field

const
  NMTMP* = 6
  NFTMP* = 4
  NBTMP* = 2

type
  TmpStruct[T,U] = object
    #mtmp*: ptr ptr QDP_ColorMatrix
    mtmp: array[NMTMP, T]
    #ftmp0*: ptr ptr QDP_ColorMatrix
    ftmp0: array[NFTMP, T]
    #ftmp*: array[NFTMP, ptr ptr QDP_ColorMatrix]
    ftmp: array[NFTMP, array[4,Shifter[T,U]]]
    #btmp0*: ptr ptr QDP_ColorMatrix
    btmp0: array[NBTMP, T]
    #btmp*: array[NBTMP, ptr ptr QDP_ColorMatrix]
    btmp: array[NBTMP, array[4,Shifter[T,U]]]
    #lat*: ptr QDP_Lattice
    #sub*: QDP_Subset
    #nc*: cint
    #nd*: cint

#const
#  PEQM* = (2 * NC * NC)
#  EQMTM* = (NC * NC * (8 * NC - 2))
#  PEQMTM* = (NC * NC * (8 * NC))

proc set_temps(t: var TmpStruct, lo: Layout) =
  #var i: cint = 0
  #while i < NMTMP:
  for i in 0..<NMTMP:
    t.mtmp[i].new(lo)
  #var i: cint = 0
  #while i < NFTMP:
  for i in 0..<NFTMP:
    t.ftmp0[i].new(lo)
    #var j: cint = 0
    #while j < t.nd:
    for j in 0..<4:
      t.ftmp[i][j] = newShifter(t.ftmp0[i], j, 1)
  #var i: cint = 0
  #while i < NBTMP:
  for i in 0..<NBTMP:
    t.btmp0[i].new(lo)
    #var j: cint = 0
    #while j < t.nd:
    for j in 0..<4:
      t.btmp[i][j] = newShifter(t.btmp0[i], j, -1)

proc free_temps(t: var TmpStruct) =
  for i in 0..<NMTMP:
    t.mtmp[i] = nil
  for i in 0..<NFTMP:
    t.ftmp0[i] = nil
    #for j in 0..<4:
    #  t.ftmp[i][j] = nil
  for i in 0..<NBTMP:
    t.btmp0[i] = nil
    #for j in 0..<4:
    #  t.btmp[i][j] = nil

proc staples(stap: auto; top0: auto; bot: auto; left: auto; right0: auto;
             mu: int; nu: int; t: TmpStruct) =
  let
    right = t.ftmp0[0]
    rightmu = t.ftmp[0][mu]
    top = t.ftmp0[1]
    topnu = t.ftmp[1][nu]
    leftbot = t.mtmp[0]
    back = t.btmp0[0]
    backnu = t.btmp[0][nu]
    lefttopnu = t.mtmp[1]
  #QDP_M_eq_M(right, right0, t.sub)
  right := right0
  #QDP_M_eq_sM(rightmu, right, QDP_neighbor[mu], QDP_forward, t.sub)
  discard rightmu ^* right
  #QDP_M_eq_M(top, top0, t.sub)
  top := top0
  #QDP_M_eq_sM(topnu, top, QDP_neighbor[nu], QDP_forward, t.sub)
  discard topnu ^* top
  #QDP_M_eq_Ma_times_M(leftbot, left, bot, t.sub)
  leftbot := left.adj * bot
  #QDP_M_eq_M_times_M(back, leftbot, rightmu, t.sub)
  back := leftbot * rightmu.field
  #QDP_M_eq_sM(backnu, back, QDP_neighbor[nu], QDP_backward, t.sub)
  discard backnu ^* back
  #QDP_M_eq_M_times_M(lefttopnu, left, topnu, t.sub)
  lefttopnu := left * topnu.field
  #QDP_discard_M(topnu)
  #QDP_M_peq_M_times_Ma(stap, lefttopnu, rightmu, t.sub)
  stap += lefttopnu * rightmu.field.adj
  #QDP_discard_M(rightmu)
  #QDP_M_peq_M(stap, backnu, t.sub)
  stap += backnu.field
  #QDP_discard_M(backnu)
  #const
  #  STAPLES_FLOPS = (3 * EQMTM + PEQMTM + PEQM)

proc staple2fb(outmuf: auto; outmub: auto; outnuf: auto; outnub: auto;
               Umu0: auto; Unu0: auto; mu: int; nu: int; t: TmpStruct) =
  let
    Unu = t.ftmp0[0]
    Unufmu = t.ftmp[0][mu]
    Umu = t.ftmp0[1]
    Umufnu = t.ftmp[1][nu]
    UmuUnufmu = t.mtmp[0]
    backmu = t.btmp0[0]
    backmubnu = t.btmp[0][nu]
    UnuUmufnu = t.mtmp[1]
    backnu = t.btmp0[1]
    backnubmu = t.btmp[1][mu]
  ##  fmu += Unu * Umufnu * Unufmu+
  ##  fmu += Unubnu+ * Umubnu * Unufmubnu
  ##  fnu += Umu * Unufmu * Umufnu+
  ##  fnu += Umubmu+ * Unubmu * Umufnubmu
  #QDP_M_eq_M(Unu, Unu0, t.sub)
  Unu := Unu0
  #QDP_M_eq_sM(Unufmu, Unu, QDP_neighbor[mu], QDP_forward, t.sub)
  discard Unufmu ^* Unu
  #QDP_M_eq_M(Umu, Umu0, t.sub)
  Umu := Umu0
  #QDP_M_eq_sM(Umufnu, Umu, QDP_neighbor[nu], QDP_forward, t.sub)
  discard Umufnu ^* Umu
  #QDP_M_eq_M_times_M(UmuUnufmu, Umu, Unufmu, t.sub)
  UmuUnufmu := Umu * Unufmu.field
  #QDP_M_eq_Ma_times_M(backmu, Unu, UmuUnufmu, t.sub)
  backmu := Unu.adj * UmuUnufmu
  #QDP_M_eq_sM(backmubnu, backmu, QDP_neighbor[nu], QDP_backward, t.sub)
  discard backmubnu ^* backmu
  #QDP_M_eq_M_times_M(UnuUmufnu, Unu, Umufnu, t.sub)
  UnuUmufnu := Unu * Umufnu.field
  #QDP_M_eq_Ma_times_M(backnu, Umu, UnuUmufnu, t.sub)
  backnu := Umu.adj * UnuUmufnu
  #QDP_M_eq_sM(backnubmu, backnu, QDP_neighbor[mu], QDP_backward, t.sub)
  discard backnubmu ^* backnu
  #QDP_M_eq_M_times_Ma(outmuf, UnuUmufnu, Unufmu, t.sub)
  outmuf := UnuUmufnu * Unufmu.field.adj
  #QDP_M_eq_M_times_Ma(outnuf, UmuUnufmu, Umufnu, t.sub)
  outnuf := UmuUnufmu * Umufnu.field.adj
  #QDP_discard_M(Umufnu)
  #QDP_discard_M(Unufmu)
  #QDP_M_eq_M(outmub, backmubnu, t.sub)
  outmub := backmubnu.field
  #QDP_discard_M(backmubnu)
  #QDP_M_eq_M(outnub, backnubmu, t.sub)
  outnub := backnubmu.field
  #QDP_discard_M(backnubmu)
  #const
  #  STAPLE2FB_FLOPS = (6 * EQMTM)

proc stapler(outmu: auto; outnu: auto; Umu0: auto; Unu0: auto; Fmu0: auto;
             Fnu0: auto; Bmu: auto; Bnu: auto; mu: int; nu: int; t: TmpStruct) =
  let
    Unu = t.ftmp0[0]
    Unufmu = t.ftmp[0][mu]
    Umu = t.ftmp0[1]
    Umufnu = t.ftmp[1][nu]
    Fnu = t.ftmp0[2]
    Fnufmu = t.ftmp[2][mu]
    Fmu = t.ftmp0[3]
    Fmufnu = t.ftmp[3][nu]
    UmuUnufmu = t.mtmp[0]
    BmuUnufmu = t.mtmp[1]
    backmu = t.btmp0[0]
    backmubnu = t.btmp[0][nu]
    UnuUmufnu = t.mtmp[2]
    BnuUmufnu = t.mtmp[3]
    backnu = t.btmp0[1]
    backnubmu = t.btmp[1][mu]
  ##  fmu += Unu * Umufnu * Fnufmu+
  ##  fmu += Unu * Fmufnu * Unufmu+
  ##  fmu += Bnu * Umufnu * Unufmu+
  ##  fmu += Unubnu+ * Umubnu * Fnufmubnu
  ##  fmu += Unubnu+ * Bmubnu * Unufmubnu
  ##  fmu += Bnubnu+ * Umubnu * Unufmubnu
  ##  fnu += Umu * Unufmu * Fmufnu+
  ##  fnu += Umu * Fnufmu * Umufnu+
  ##  fnu += Bmu * Unufmu * Umufnu+
  ##  fnu += Umubmu+ * Unubmu * Fmufnubmu
  ##  fnu += Umubmu+ * Bnubmu * Umufnubmu
  ##  fnu += Bmubmu+ * Unubmu * Umufnubmu
  #QDP_M_eq_M(Unu, Unu0, t.sub)
  Unu := Unu0
  #QDP_M_eq_sM(Unufmu, Unu, QDP_neighbor[mu], QDP_forward, t.sub)
  discard Unufmu ^* Unu
  #QDP_M_eq_M(Umu, Umu0, t.sub)
  Umu := Umu0
  #QDP_M_eq_sM(Umufnu, Umu, QDP_neighbor[nu], QDP_forward, t.sub)
  discard Umufnu ^* Umu
  #QDP_M_eq_M(Fnu, Fnu0, t.sub)
  Fnu := Fnu0
  #QDP_M_eq_sM(Fnufmu, Fnu, QDP_neighbor[mu], QDP_forward, t.sub)
  discard Fnufmu ^* Fnu
  #QDP_M_eq_M(Fmu, Fmu0, t.sub)
  Fmu := Fmu0
  #QDP_M_eq_sM(Fmufnu, Fmu, QDP_neighbor[nu], QDP_forward, t.sub)
  discard Fmufnu ^* Fmu
  #QDP_M_eq_M_times_M(UmuUnufmu, Umu, Unufmu, t.sub)
  UmuUnufmu := Umu * Unufmu.field
  #QDP_M_eq_Ma_times_M(backmu, Bnu, UmuUnufmu, t.sub)
  backmu := Bnu.adj * UmuUnufmu
  #QDP_M_eq_M_times_M(BmuUnufmu, Bmu, Unufmu, t.sub)
  BmuUnufmu := Bmu * Unufmu.field
  #QDP_M_peq_M_times_M(BmuUnufmu, Umu, Fnufmu, t.sub)
  BmuUnufmu += Umu * Fnufmu.field
  #QDP_M_peq_Ma_times_M(backmu, Unu, BmuUnufmu, t.sub)
  backmu += Unu.adj * BmuUnufmu
  #QDP_M_eq_sM(backmubnu, backmu, QDP_neighbor[nu], QDP_backward, t.sub)
  discard backmubnu ^* backmu
  #QDP_M_eq_M_times_M(UnuUmufnu, Unu, Umufnu, t.sub)
  UnuUmufnu := Unu * Umufnu.field
  #QDP_M_eq_Ma_times_M(backnu, Bmu, UnuUmufnu, t.sub)
  backnu := Bmu.adj * UnuUmufnu
  #QDP_M_eq_M_times_M(BnuUmufnu, Bnu, Umufnu, t.sub)
  BnuUmufnu := Bnu * Umufnu.field
  #QDP_M_peq_M_times_M(BnuUmufnu, Unu, Fmufnu, t.sub)
  BnuUmufnu += Unu * Fmufnu.field
  #QDP_M_peq_Ma_times_M(backnu, Umu, BnuUmufnu, t.sub)
  backnu += Umu.adj * BnuUmufnu
  #QDP_M_eq_sM(backnubmu, backnu, QDP_neighbor[mu], QDP_backward, t.sub)
  discard backnubmu ^* backnu
  #QDP_M_peq_M_times_Ma(outmu, UnuUmufnu, Fnufmu, t.sub)
  outmu += UnuUmufnu * Fnufmu.field.adj
  #QDP_M_peq_M_times_Ma(outmu, BnuUmufnu, Unufmu, t.sub)
  outmu += BnuUmufnu * Unufmu.field.adj
  #QDP_M_peq_M(outmu, backmubnu, t.sub)
  outmu += backmubnu.field
  #QDP_M_peq_M_times_Ma(outnu, UmuUnufmu, Fmufnu, t.sub)
  outnu += UmuUnufmu * Fmufnu.field.adj
  #QDP_M_peq_M_times_Ma(outnu, BmuUnufmu, Umufnu, t.sub)
  outnu += BmuUnufmu * Umufnu.field.adj
  #QDP_M_peq_M(outnu, backnubmu, t.sub)
  outnu += backnubmu.field
  #QDP_discard_M(Umufnu)
  #QDP_discard_M(Unufmu)
  #QDP_discard_M(Fmufnu)
  #QDP_discard_M(Fnufmu)
  #QDP_discard_M(backmubnu)
  #QDP_discard_M(backnubmu)
  #const
  #  STAPLER_FLOPS = (6 * EQMTM + 8 * PEQMTM + 2 * PEQM)

proc staplep*(outmu: auto; outnu: auto; Umu0: auto; Unu0: auto; Xmu0: auto; Xnu0: auto;
              nx: int; mu: int; nu: int; t: TmpStruct) =
  let
    XXnumu = t.mtmp[0]
    XXmunu = t.mtmp[1]
    XXdnumu = t.mtmp[2]
    XXdmunu = t.mtmp[3]
    XdXnumu = t.mtmp[4]
    XdXmunu = t.mtmp[5]
    Xnu = t.ftmp0[0]
    Xnufmu = t.ftmp[0][mu]
    Xmu = t.ftmp0[1]
    Xmufnu = t.ftmp[1][nu]
    Unu = t.ftmp0[2]
    Unufmu = t.ftmp[2][mu]
    Umu = t.ftmp0[3]
    Umufnu = t.ftmp[3][nu]
    backmu = t.btmp0[0]
    backmubnu = t.btmp[0][nu]
    backnu = t.btmp0[1]
    backnubmu = t.btmp[1][mu]
  ##  fmu += Unu * Xmufnu * Xnufmu+
  ##  fmu += Unubnu+ * Xmubnu * Xnufmubnu
  ##  fmu += Xnu * Xmufnu * Unufmu+
  ##  fmu += Xnubnu+ * Xmubnu * Unufmubnu
  ##  fnu += Umu * Xnufmu * Xmufnu+
  ##  fnu += Umubmu+ * Xnubmu * Xmufnubmu
  ##  fnu += Xmu * Xnufmu * Umufnu+
  ##  fnu += Xmubmu+ * Xnubmu * Umufnubmu
  #QDP_M_eq_zero(XXnumu, t.sub)
  XXnumu := 0
  #QDP_M_eq_zero(XXmunu, t.sub)
  XXmunu := 0
  #QDP_M_eq_zero(XXdnumu, t.sub)
  XXdnumu := 0
  #QDP_M_eq_zero(XXdmunu, t.sub)
  XXdmunu := 0
  #QDP_M_eq_zero(XdXnumu, t.sub)
  XdXnumu := 0
  #QDP_M_eq_zero(XdXmunu, t.sub)
  XdXmunu := 0
  #var i: cint = 0
  #while i < nx:
  for i in 0..<nx:
    #QDP_M_eq_M(Xnu, Xnu0[i], t.sub)
    Xnu := Xnu0[i]
    #QDP_M_eq_sM(Xnufmu, Xnu, QDP_neighbor[mu], QDP_forward, t.sub)
    discard Xnufmu ^* Xnu
    #QDP_M_eq_M(Xmu, Xmu0[i], t.sub)
    Xmu := Xmu0[i]
    #QDP_M_eq_sM(Xmufnu, Xmu, QDP_neighbor[nu], QDP_forward, t.sub)
    discard Xmufnu ^* Xmu
    #QDP_M_peq_Ma_times_M(XdXnumu, Xnu, Xmu, t.sub)
    XdXnumu += Xnu.adj * Xmu
    #QDP_M_peq_Ma_times_M(XdXmunu, Xmu, Xnu, t.sub)
    XdXmunu += Xmu.adj * Xnu
    #QDP_M_peq_M_times_M(XXnumu, Xnu, Xmufnu, t.sub)
    XXnumu += Xnu * Xmufnu.field
    #QDP_M_peq_M_times_M(XXmunu, Xmu, Xnufmu, t.sub)
    XXmunu += Xmu * Xnufmu.field
    #QDP_M_peq_M_times_Ma(XXdnumu, Xnufmu, Xmufnu, t.sub)
    XXdnumu += Xnufmu.field * Xmufnu.field.adj
    #QDP_M_peq_M_times_Ma(XXdmunu, Xmufnu, Xnufmu, t.sub)
    XXdmunu += Xmufnu.field * Xnufmu.field.adj
    #QDP_discard_M(Xmufnu)
    #QDP_discard_M(Xnufmu)
  #QDP_M_eq_M(Unu, Unu0, t.sub)
  Unu := Unu0
  #QDP_M_eq_sM(Unufmu, Unu, QDP_neighbor[mu], QDP_forward, t.sub)
  discard Unufmu ^* Unu
  #QDP_M_eq_M(Umu, Umu0, t.sub)
  Umu := Umu0
  #QDP_M_eq_sM(Umufnu, Umu, QDP_neighbor[nu], QDP_forward, t.sub)
  discard Umufnu ^* Umu
  #QDP_M_eq_Ma_times_M(backmu, Unu, XXmunu, t.sub)
  backmu := Unu.adj * XXmunu
  #QDP_M_peq_M_times_M(backmu, XdXnumu, Unufmu, t.sub)
  backmu += XdXnumu * Unufmu.field
  #QDP_M_eq_sM(backmubnu, backmu, QDP_neighbor[nu], QDP_backward, t.sub)
  discard backmubnu ^* backmu
  #QDP_M_eq_Ma_times_M(backnu, Umu, XXnumu, t.sub)
  backnu := Umu.adj * XXnumu
  #QDP_M_peq_M_times_M(backnu, XdXmunu, Umufnu, t.sub)
  backnu += XdXmunu * Umufnu.field
  #QDP_M_eq_sM(backnubmu, backnu, QDP_neighbor[mu], QDP_backward, t.sub)
  discard backnubmu ^* backnu
  #QDP_M_peq_M_times_M(outmu, Unu, XXdmunu, t.sub)
  outmu += Unu * XXdmunu
  #QDP_M_peq_M_times_Ma(outmu, XXnumu, Unufmu, t.sub)
  outmu += XXnumu * Unufmu.field.adj
  #QDP_M_peq_M(outmu, backmubnu, t.sub)
  outmu += backmubnu.field
  #QDP_M_peq_M_times_M(outnu, Umu, XXdnumu, t.sub)
  outnu += Umu * XXdnumu
  #QDP_M_peq_M_times_Ma(outnu, XXmunu, Umufnu, t.sub)
  outnu += XXmunu * Umufnu.field.adj
  #QDP_M_peq_M(outnu, backnubmu, t.sub)
  outnu += backnubmu.field
  #QDP_discard_M(Umufnu)
  #QDP_discard_M(Unufmu)
  #QDP_discard_M(backmubnu)
  #QDP_discard_M(backnubmu)
  #template STAPLEP_FLOPS(n: untyped): untyped =
  #  (2 * EQMTM + (6 * n + 6) * PEQMTM + 2 * PEQM)

proc staplep2*(outmu: auto; Unu0: auto; Xmu0: auto; Xnu0: auto;
               nx: int; mu: int; nu: int; t: TmpStruct) =
  let
    XXnumu = t.mtmp[0]
    XXmunu = t.mtmp[1]
    XXdmunu = t.mtmp[2]
    XdXnumu = t.mtmp[3]
    Xnu = t.ftmp0[0]
    Xnufmu = t.ftmp[0][mu]
    Xmu = t.ftmp0[1]
    Xmufnu = t.ftmp[1][nu]
    Unu = t.ftmp0[2]
    Unufmu = t.ftmp[2][mu]
    backmu = t.btmp0[0]
    backmubnu = t.btmp[0][nu]
  ##  fmu += Unu * Xmufnu * Xnufmu+
  ##  fmu += Unubnu+ * Xmubnu * Xnufmubnu
  ##  fmu += Xnu * Xmufnu * Unufmu+
  ##  fmu += Xnubnu+ * Xmubnu * Unufmubnu
  #var i: cint = 0
  #while i < nx:
  for i in 0..<nx:
    #QDP_M_eq_M(Xnu, Xnu0[i], t.sub)
    Xnu := Xnu0[i]
    #QDP_M_eq_sM(Xnufmu, Xnu, QDP_neighbor[mu], QDP_forward, t.sub)
    discard Xnufmu ^* Xnu
    #QDP_M_eq_M(Xmu, Xmu0[i], t.sub)
    Xmu := Xmu0[i]
    #QDP_M_eq_sM(Xmufnu, Xmu, QDP_neighbor[nu], QDP_forward, t.sub)
    discard Xmufnu ^* Xmu
    if i == 0:
      #QDP_M_eq_Ma_times_M(XdXnumu, Xnu, Xmu, t.sub)
      XdXnumu := Xnu.adj * Xmu
      #QDP_M_eq_M_times_M(XXnumu, Xnu, Xmufnu, t.sub)
      XXnumu := Xnu * Xmufnu.field
      #QDP_M_eq_M_times_M(XXmunu, Xmu, Xnufmu, t.sub)
      XXmunu := Xmu * Xnufmu.field
      #QDP_M_eq_M_times_Ma(XXdmunu, Xmufnu, Xnufmu, t.sub)
      XXdmunu := Xmufnu.field * Xnufmu.field.adj
    else:
      #QDP_M_peq_Ma_times_M(XdXnumu, Xnu, Xmu, t.sub)
      XdXnumu += Xnu.adj * Xmu
      #QDP_M_peq_M_times_M(XXnumu, Xnu, Xmufnu, t.sub)
      XXnumu += Xnu * Xmufnu.field
      #QDP_M_peq_M_times_M(XXmunu, Xmu, Xnufmu, t.sub)
      XXmunu += Xmu * Xnufmu.field
      #QDP_M_peq_M_times_Ma(XXdmunu, Xmufnu, Xnufmu, t.sub)
      XXdmunu += Xmufnu.field * Xnufmu.field.adj
    #QDP_discard_M(Xmufnu)
    #QDP_discard_M(Xnufmu)
    #inc(i)
  #QDP_M_eq_M(Unu, Unu0, t.sub)
  Unu := Unu0
  #QDP_M_eq_sM(Unufmu, Unu, QDP_neighbor[mu], QDP_forward, t.sub)
  discard Unufmu ^* Unu
  #QDP_M_eq_Ma_times_M(backmu, Unu, XXmunu, t.sub)
  backmu := Unu.adj * XXmunu
  #QDP_M_peq_M_times_M(backmu, XdXnumu, Unufmu, t.sub)
  backmu += XdXnumu * Unufmu.field
  #QDP_M_eq_sM(backmubnu, backmu, QDP_neighbor[nu], QDP_backward, t.sub)
  discard backmubnu ^* backmu
  #QDP_M_peq_M_times_M(outmu, Unu, XXdmunu, t.sub)
  outmu += Unu * XXdmunu
  #QDP_M_peq_M_times_Ma(outmu, XXnumu, Unufmu, t.sub)
  outmu += XXnumu * Unufmu.field.adj
  #QDP_M_peq_M(outmu, backmubnu, t.sub)
  outmu += backmubnu.field
  #QDP_discard_M(Unufmu)
  #QDP_discard_M(backmubnu)
  #template STAPLEP2_FLOPS(n: untyped): untyped =
  #  (5 * EQMTM + (4 * (n - 1) + 3) * PEQMTM + PEQM)

proc symanzik1loopDeriv*(coeffs: GaugeActionCoeffs;
                         links: auto; deriv: auto;
                         eps: float; info: var PerfInfo;) =
  tic("symanzik1loopDeriv")
  template G(x: untyped): untyped = links[x]
  type
    GF = type links[0]
    GM = type links[0][0]
    Shft = Shifter[GF,GM]
  const
    nd = 4
  let
    nc = links[0][0].getNc
    lo = links[0].l
    #nd = lat.nDim
    fac = - eps / nc
    plaq = fac * coeffs.plaq
    rect = fac * coeffs.rect
    pgm = fac * coeffs.pgm
    adpl = 2 * (fac/nc) * coeffs.adjplaq
  var
    #dtime: cdouble = QOP_time()
    nflops = newThreadSingle(0)
    t: TmpStruct[GF,GM]
  #t.nc = QOP_Nc
  #t.nd = nd
  #t.sub = sub
  #t.lat = lat
  set_temps(t, lo)
  toc("set_temps")
  var
    tforce: array[nd, GF]
    stplf: array[nd, array[nd, GF]]
    stplb: array[nd, array[nd, GF]]
    tmat: array[nd, GF]
    tc = lo.Complex
  #var mu: cint = 0
  #while mu < nd:
  for mu in 0..<nd:
    tforce[mu].new(lo)
  threads:
    for mu in 0..<nd:
      tforce[mu] := 0
  if rect != 0.0 or pgm != 0.0:
    #var mu: cint = 0
    #while mu < nd:
    for mu in 0..<nd:
      tmat[mu].new(lo)
  #var mu: cint = 1
  #while mu < nd:
  for mu in 0..<nd:
    #var nu: cint = 0
    #while nu < mu:
    for nu in 0..<mu:
      stplf[mu][nu].new(lo)
      stplb[mu][nu].new(lo)
      stplf[nu][mu].new(lo)
      stplb[nu][mu].new(lo)
  toc("setup")
  threads:
    for mu in 0..<nd:
      #var nu: cint = 0
      #while nu < mu:
      for nu in 0..<mu:
        staple2fb(stplf[mu][nu], stplb[mu][nu], stplf[nu][mu], stplb[nu][mu],
                  G(mu), G(nu), mu, nu, t)
        #nflops += STAPLE2FB_FLOPS
        if adpl != 0.0:
          #var z: QLA_Complex
          #QLA_c_eq_r(z, plaq div adpl)
          #QDP_C_eq_c(tc, addr(z), sub)
          #QDP_C_peq_M_dot_M(tc, stplf[mu][nu], G(mu), sub)
          #QDP_C_eq_r_times_C(tc, addr(adpl), tc, sub)
          for e in lo:
            tc[e] := plaq + adpl*dot(stplf[mu][nu][e], G(mu)[e])
          #QDP_M_peq_C_times_M(tforce[mu], tc, stplf[mu][nu], sub)
          tforce[mu] += tc * stplf[mu][nu]
          #QDP_C_eq_c(tc, addr(z), sub)
          #QDP_C_peq_M_dot_M(tc, stplb[mu][nu], G(mu), sub)
          #QDP_C_eq_r_times_C(tc, addr(adpl), tc, sub)
          for e in lo:
            tc[e] := plaq + adpl*dot(stplb[mu][nu][e], G(mu)[e])
          #QDP_M_peq_C_times_M(tforce[mu], tc, stplb[mu][nu], sub)
          tforce[mu] += tc * stplb[mu][nu]
          #QDP_C_eq_c(tc, addr(z), sub)
          #QDP_C_peq_M_dot_M(tc, stplf[nu][mu], G(nu), sub)
          #QDP_C_eq_r_times_C(tc, addr(adpl), tc, sub)
          for e in lo:
            tc[e] := plaq + adpl*dot(stplf[nu][mu][e], G(nu)[e])
          #QDP_M_peq_C_times_M(tforce[nu], tc, stplf[nu][mu], sub)
          tforce[nu] += tc * stplf[nu][mu]
          #QDP_C_eq_c(tc, addr(z), sub)
          #QDP_C_peq_M_dot_M(tc, stplb[nu][mu], G(nu), sub)
          #QDP_C_eq_r_times_C(tc, addr(adpl), tc, sub)
          for e in lo:
            tc[e] := plaq + adpl*dot(stplb[nu][mu][e], G(nu)[e])
          #QDP_M_peq_C_times_M(tforce[nu], tc, stplb[nu][mu], sub)
          tforce[nu] += tc * stplb[nu][mu]
          #QDP_destroy_C(tc)
          #nflops += 4 * (16 * NC * NC + 2)
        elif plaq != 0.0:
          #QDP_M_peq_r_times_M(tforce[mu], addr(plaq), stplf[mu][nu], sub)
          tforce[mu] += plaq * stplf[mu][nu]
          #QDP_M_peq_r_times_M(tforce[mu], addr(plaq), stplb[mu][nu], sub)
          tforce[mu] += plaq * stplb[mu][nu]
          #QDP_M_peq_r_times_M(tforce[nu], addr(plaq), stplf[nu][mu], sub)
          tforce[nu] += plaq * stplf[nu][mu]
          #QDP_M_peq_r_times_M(tforce[nu], addr(plaq), stplb[nu][mu], sub)
          tforce[nu] += plaq * stplb[nu][mu]
          #nflops += 4 * 4 * NC * NC
    toc("plaq & adjplaq")
    if rect != 0.0:
      #var mu: cint = 0
      #while mu < nd:
      for mu in 0..<nd:
        #QDP_M_eq_zero(tmat[mu], sub)
        tmat[mu] := 0
        #inc(mu)
      #var mu: cint = 1
      #while mu < nd:
      for mu in 1..<nd:
        #var nu: cint = 0
        #while nu < mu:
        for nu in 0..<mu:
          stapler(tmat[mu], tmat[nu], G(mu), G(nu), stplf[mu][nu], stplf[nu][mu],
                  stplb[mu][nu], stplb[nu][mu], mu, nu, t)
          #nflops += STAPLER_FLOPS
          #inc(nu)
        #inc(mu)
      #var mu: cint = 0
      #while mu < nd:
      for mu in 0..<nd:
        #QDP_M_peq_r_times_M(tforce[mu], addr(rect), tmat[mu], sub)
        tforce[mu] += rect * tmat[mu]
        #inc(mu)
      #nflops += nd * 4 * NC * NC
    toc("rect")
    if pgm != 0.0:
      if nd != 4:
        #var mu: cint = 0
        #while mu < nd:
        for mu in 0..<nd:
          #QDP_M_eq_zero(tmat[0], sub)
          tmat[0] := 0
          #var nu: cint = 0
          #while nu < nd:
          for nu in 0..<nd:
            if nu == mu: continue
            #var rho = nu + 1
            #while rho < nd:
            for rho in (nu+1)..<nd:
              if rho == mu or rho == nu: continue
              staples(tmat[0], stplf[mu][rho], stplf[mu][rho], stplf[nu][rho], G(nu),
                      mu, nu, t)
              staples(tmat[0], stplf[mu][rho], stplf[mu][rho], G(nu), stplf[nu][rho],
                      mu, nu, t)
              staples(tmat[0], stplb[mu][rho], stplb[mu][rho], stplb[nu][rho], G(nu),
                      mu, nu, t)
              staples(tmat[0], stplb[mu][rho], stplb[mu][rho], G(nu), stplb[nu][rho],
                      mu, nu, t)
              #nflops += 4 * STAPLES_FLOPS
          #QDP_M_peq_r_times_M(tforce[mu], addr(pgm), tmat[0], sub)
          tforce[mu] += pgm * tmat[0]
          #nflops += 4 * NC * NC
      else:
        ##  nd == 4
        #var mu: cint = 0
        #while mu < 4:
        for mu in 0..<4:
          #QDP_M_eq_zero(tmat[mu], sub)
          tmat[mu] := 0
        var
          mu: int
          nu: int
          rho: int
        var
          Xmu: array[4, GF]
          Xnu: array[4, GF]
        mu = 0
        nu = 1
        Xmu[0] = stplf[mu][2]
        Xmu[1] = stplf[mu][3]
        Xmu[2] = stplb[mu][2]
        Xmu[3] = stplb[mu][3]
        Xnu[0] = stplf[nu][2]
        Xnu[1] = stplf[nu][3]
        Xnu[2] = stplb[nu][2]
        Xnu[3] = stplb[nu][3]
        staplep(tmat[mu], tmat[nu], G(mu), G(nu), Xmu, Xnu, 4, mu, nu, t)
        #nflops += STAPLEP_FLOPS(4)
        mu = 2
        nu = 3
        Xmu[0] = stplf[mu][0]
        Xmu[1] = stplf[mu][1]
        Xmu[2] = stplb[mu][0]
        Xmu[3] = stplb[mu][1]
        Xnu[0] = stplf[nu][0]
        Xnu[1] = stplf[nu][1]
        Xnu[2] = stplb[nu][0]
        Xnu[3] = stplb[nu][1]
        staplep(tmat[mu], tmat[nu], G(mu), G(nu), Xmu, Xnu, 4, mu, nu, t)
        #nflops += STAPLEP_FLOPS(4)
        mu = 0
        nu = 2
        rho = 3
        Xmu[0] = stplf[mu][rho]
        Xmu[1] = stplb[mu][rho]
        Xnu[0] = stplf[nu][rho]
        Xnu[1] = stplb[nu][rho]
        staplep2(tmat[mu], G(nu), Xmu, Xnu, 2, mu, nu, t)
        mu = 1
        nu = 2
        rho = 3
        Xmu[0] = stplf[mu][rho]
        Xmu[1] = stplb[mu][rho]
        Xnu[0] = stplf[nu][rho]
        Xnu[1] = stplb[nu][rho]
        staplep2(tmat[mu], G(nu), Xmu, Xnu, 2, mu, nu, t)
        mu = 2
        nu = 0
        rho = 1
        Xmu[0] = stplf[mu][rho]
        Xmu[1] = stplb[mu][rho]
        Xnu[0] = stplf[nu][rho]
        Xnu[1] = stplb[nu][rho]
        staplep2(tmat[mu], G(nu), Xmu, Xnu, 2, mu, nu, t)
        mu = 3
        nu = 0
        rho = 1
        Xmu[0] = stplf[mu][rho]
        Xmu[1] = stplb[mu][rho]
        Xnu[0] = stplf[nu][rho]
        Xnu[1] = stplb[nu][rho]
        staplep2(tmat[mu], G(nu), Xmu, Xnu, 2, mu, nu, t)
        #nflops += 4 * STAPLEP2_FLOPS(2)
        #var mu: cint = 0
        #while mu < 4:
        for mu in 0..<4:
          #QDP_M_peq_r_times_M(tforce[mu], addr(pgm), tmat[mu], sub)
          tforce[mu] += pgm * tmat[mu]
        #nflops += 4 * 4 * NC * NC
      ##  nd!=4
    toc("pgm")
    #var mu: cint = 0
    #while mu < nd:
    for mu in 0..<nd:
      #QDP_M_peq_M(deriv[mu], tforce[mu], sub)
      deriv[mu] += tforce[mu]
      #inc(mu)
    #nflops += nd * PEQM
  #if rect or pgm:
  #  var mu: cint = 0
  #  while mu < nd:
  #    QDP_destroy_M(tmat[mu])
  #    inc(mu)
  #var mu: cint = 0
  #while mu < nd:
  #  QDP_destroy_M(tforce[mu])
  #  var nu: cint = 0
  #  while nu < nd:
  #    if nu == mu:
  #      inc(nu)
  #      inc(mu)
  #      continue
  #    QDP_destroy_M(stplf[mu][nu])
  #    QDP_destroy_M(stplb[mu][nu])
  #    inc(nu)
  #  inc(mu)
  free_temps(t)
  ## double nflop = 96720 - 4*(24+18);
  #info.final_sec = QOP_time() - dtime
  #info.final_flop = nflops * QDP_sites_on_node
  #info.status = QOP_SUCCESS

#[
proc QOP_symanzik_1loop_gauge_force_qdp*(info: ptr QOP_info_t;
                                        links: ptr ptr QDP_ColorMatrix;
                                        force: ptr ptr QDP_ColorMatrix;
                                        coeffs: ptr QOP_gauge_coeffs_t;
                                        eps: QLA_Real) =
  const
    NC = QDP_get_nc(links[0])
  var dtime: cdouble = QOP_time()
  var lat: ptr QDP_Lattice = QDP_get_lattice_M(links[0])
  var nd: cint = QDP_ndim_L(lat)
  var sub: QDP_Subset = QDP_all_L(lat)
  var deriv: array[nd, ptr QDP_ColorMatrix]
  var mu: cint = 0
  while mu < nd:
    deriv[mu] = QDP_create_M_L(lat)
    QDP_M_eq_zero(deriv[mu], sub)
    inc(mu)
  QOP_symanzik_1loop_gauge_deriv_qdp(info, links, deriv, coeffs, eps, 0)
  var mtmp: ptr QDP_ColorMatrix = QDP_create_M_L(lat)
  when defined(CHKSUM):
    QLA_ColorMatrix(qcm)
    var
      det: QLA_Complex
      chk: QLA_Complex
    QLA_c_eq_r(chk, 0)
  var trace: QLA_Real = 0
  var mu: cint = 0
  while mu < nd:
    QDP_M_eq_M_times_Ma(mtmp, G(mu), deriv[mu], sub)
    if QOP_common.verbosity == QOP_VERB_DEBUG:
      QLA_ColorMatrix(tm)
      var tr: QLA_Real
      QDP_m_eq_sum_M(addr(tm), mtmp, sub)
      QLA_R_eq_re_trace_M(addr(tr), addr(tm))
      inc(trace, tr)
    QDP_M_eq_antiherm_M(deriv[mu], mtmp, sub)
    QDP_M_peq_M(force[mu], deriv[mu], sub)
    when defined(CHKSUM):
      QDP_m_eq_sum_M(addr(qcm), force.force[mu], sub)
      QLA_C_eq_det_M(addr(det), addr(qcm))
      QLA_c_peq_c(chk, det)
    inc(mu)
  when defined(CHKSUM):
    QOP_printf0("chksum: %g %g\n", QLA_real(chk), QLA_imag(chk))
  if QOP_common.verbosity == QOP_VERB_DEBUG:
    QOP_printf0("re trace: %g\n", trace)
  inc(info.final_flop, nd * (EQMTM + 2 * NC * NC + PEQM) * QDP_subset_len(sub))
  QDP_destroy_M(mtmp)
  var mu: cint = 0
  while mu < nd:
    QDP_destroy_M(deriv[mu])
    inc(mu)
  info.final_sec = QOP_time() - dtime

proc QOP_symanzik_1loop_gauge_force*(info: ptr QOP_info_t;
                                    gauge: ptr QOP_GaugeField;
                                    force: ptr QOP_Force;
                                    coeffs: ptr QOP_gauge_coeffs_t; eps: REAL) =
  QOP_symanzik_1loop_gauge_force_qdp(info, gauge.links, force.force, coeffs, eps)
]#

when isMainModule:
  import qex, physics/qcdTypes
  import strformat
  qexInit()
  #var defaultGaugeFile = "l88.scidac"
  let defaultLat = @[8,8,8,8]
  #defaultSetup()

  var
    (lo, g, r) = setupLattice(defaultLat)
    fl = lo.newGauge()
    fl2 = lo.newGauge()
    ch = lo.newGauge()
    gf = lo.newGauge()
    gb = lo.newGauge()
    dg = lo.newGauge()
    dg2 = lo.newGauge()
    fd = lo.newGauge()
    ld = lo.newGauge()
    gc = GaugeActionCoeffs(plaq:1.0)
    info: PerfInfo
    eps = floatParam("eps", 1e-5)

  #g.unit
  g.random r
  #dg.gaussian r
  dg2.randomTAH r
  var dgn2 = 0.0
  for mu in 0..<g.len:
    for e in g[mu]:
      #dg[mu][e] *= eps
      #gf[mu][e] := g[mu][e] + dg[mu][e]
      #gb[mu][e] := g[mu][e] - dg[mu][e]
      dg2[mu][e] *= eps
      gf[mu][e] := exp(dg2[mu][e]) * g[mu][e]
      gb[mu][e] := exp(-dg2[mu][e]) * g[mu][e]
      dg[mu][e] := gf[mu][e] - gb[mu][e]
    dgn2 += dg[mu].norm2
  echo "dgn2: ", dgn2

  var nfail = 0
  proc check(da: float, tol: float) =
    var ds = 0.0
    for mu in 0..3:
      ds += redot(dg[mu], fd[mu])
    let r = (da-ds)/da
    #echo "Checking ", name
    echo &"  da {da}  ds {ds}  rel {r}"
    let rlim = 10*tol*dgn2*dgn2
    if abs(r) > rlim:
      inc nfail
      echo &"> ERROR rel error |{r}| > {rlim}"

  proc checkS(tol: float) =
    echo "Checking ", gc
    for mu in 0..<fd.len:
      fd[mu] := 0
    var a = gc.symanzik1loopAction(g, info)
    info.clear
    resetTimers()
    a = gc.symanzik1loopAction(gb, info)
    echo "  ", info
    let a2 = gc.symanzik1loopAction(gf, info)
    echo "  ", info
    gc.symanzik1loopDeriv(g, fd, 1, info)
    check(a2.space+a2.time-a.space-a.time, tol)

  template docheck(pl,rc,pg,ad:float) =
    gc = GaugeActionCoeffs(plaq:pl,rect:rc,pgm:pg,adjplaq:ad)
    checkS(1)

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

  if nfail > 0:
    echo "*** number of failures: ", nfail
    qexError("failed")
  else:
    echo "all tests passed"

  echoProf()
  qexFinalize()
