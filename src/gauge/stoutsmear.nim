import base
import layout
import gauge

type StoutSmear*[G] = object
  alpha*: float
  gf*,f*,expaf*,ds*,cg*:G
  # add other temporaries

proc newStoutSmear*(l:Layout, alpha:float):auto =
  # we only save a reference in gf, so no allocation to gf here.
  type G = type(l.newGauge)
  StoutSmear[G](alpha:alpha, f:l.newGauge, expaf:l.newGauge, ds:l.newGauge, cg:l.newGauge)

proc smear*[G](ss:var StoutSmear, gf:G, fl:G) =
  let
    alpha = -ss.alpha  # negative from gaugeForce
    f = ss.f
    expaf = ss.expaf
    ds = ss.ds
  ss.gf = gf

  gaugeActionDeriv(GaugeActionCoeffs(plaq:1.0), gf, ds)
  threads:
    for mu in 0..<f.len:
      for e in f[mu]:
        let s = gf[mu][e]*ds[mu][e].adj
        var t{.noinit.}: evalType(f[mu][e])
        t.projectTAH s
        f[mu][e] := t
        t := exp(alpha*t)
        expaf[mu][e] := t
        fl[mu][e] := t*gf[mu][e]

# backpropagation
# z = f(g(h(x)))
# dz/dx_i = df/dg_j dg_j/dh_k dh_k/dx_i
# Dg_j = df/dg_j, Dh_k = Dg_j dg_j/dh_k, Dx_i = Dh_k dh_k/dx_i
# convention: z = 1/2 Tr(g s† + s g†), dz/dg_ij = s†_ji = D[g]†

proc gaugeForceDeriv*[G](gf:G, deriv:G, chain:G, f:G, cg:G) =
  ## gf: gauge field in
  ## deriv: derivative out
  ## chain: back propagated derivative in
  ## f: action derivative in
  ## cg: scratch space for intermediate field
  # f_x = T_a ReTr(T_a g_x sum(staples))
  # deriv_x ← d/d(g_x)[f_x] chain_x† + sum_y d/d(g_x)[f_y] chain_y†
  let nd = gf[0].l.nDim

  # d/d(g_x)[f_x] c_x† = d/d(g_x)[T_a ReTr(T_a g_x B)] c_x†
  #                    = 1/2 d/d(g_x)_ij[(g_x B T_a + c.c.)] T_akl c_x†_lk
  #                    = 1/2 (B T_a)†_ji Tr(T_a c_x†)
  # redefine the chain, add the complex conjugate back
  # D[g_x]† = (B projTAH(c_x†))†, B: x+ν→x
  # D[g_x] = projTAH(c_x) B†
  threads:
    for mu in 0..<nd:
      for e in gf[mu]:
        var t{.noinit.}: evalType(chain[mu][e])
        t.projectTAH chain[mu][e]
        deriv[mu][e] := t*f[mu][e]
        cg[mu][e] := t.adj*gf[mu][e]

  # d/d(g_x)[f_y] c_y† = d/d(g_x)[T_a ReTr(T_a g_y A g_x B)] c_y†
  #                    = 1/2 d/d(g_x)_ij[(g_x B T_a g_y A + c.c.)] T_akl c_y†_lk
  #                    = 1/2 (B T_a g_y A)†_ji Tr(T_a c_y†)
  # redefine the chain, add the complex conjugate back
  # D[g_x]† = (B projTAH(c_y†) g_y A)†, A: y+μ→x, B: x+ν→y
  # D[g_x] = A† (projTAH(c_y†) g_y)†  B†, staples: x → y+μ →[g_y† cTAH_μ] y → x+ν, with cTAH_μ insertion
  const nc = gf[0][0].nrows
  let cp = 1.0 / float(nc)    # normalization follows gaugeForce
  let
    t = newTransporters(gf, gf[0], 1)
    td = newTransporters(gf, gf[0], -1)
    ct = newTransporters(cg, cg[0], 1)
    ctd = newTransporters(cg, cg[0], -1)
  threads:
    for mu in 0..<nd:
      for nu in 0..<nd:
        if nu==mu: continue
        discard t[nu] ^* gf[mu]
        shiftExpr(t[mu].sb, deriv[mu][ir] += cp * t[nu].field[ir]*adj(it), cg[nu][ix])
        discard t[nu] ^* cg[mu]
        shiftExpr(t[mu].sb, deriv[mu][ir] += cp * t[nu].field[ir]*adj(it), gf[nu][ix])
        discard ct[nu] ^* gf[mu]
        shiftExpr(t[mu].sb, deriv[mu][ir] += cp * ct[nu].field[ir]*adj(it), gf[nu][ix])
        deriv[mu] += cp * td[nu] ^* t[mu] ^* cg[nu]
        deriv[mu] += cp * td[nu] ^* ct[mu] ^* gf[nu]
        deriv[mu] += cp * ctd[nu] ^* t[mu] ^* gf[nu]

proc smearDeriv*[G](ss:var StoutSmear, deriv:G, chain:G) =
  ## gf: gauge field in
  ## deriv: derivative out
  ## chain: back propagated derivative in
  ## chain and deriv much be distinct fields
  let
    alpha = -ss.alpha  # negative from gaugeForce
    gf = ss.gf
    f = ss.f
    expaf = ss.expaf

  # convention:
  # D[exp(a f)]† = d/d(exp(a f))[exp(a f) g] c†
  # D[exp(a f)]† = g c†
  # D[a f]† = d/d(a f)[exp(a f)] D[exp(a f)]†
  # D[f]† = d/d(f)[a f] D[a f]†
  threads:
    for mu in 0..<deriv.len:
      for e in deriv[mu]:
        deriv[mu][e] := alpha*expDeriv(alpha*f[mu][e], chain[mu][e]*gf[mu][e].adj)
  gaugeForceDeriv(ss.gf, deriv, deriv, ss.ds, ss.cg)

  # D[g]† = d/d(g)[exp(a f) g] c† = c† exp(a f)
  threads:
    for mu in 0..<deriv.len:
      for e in deriv[mu]:
        deriv[mu][e] += expaf[mu][e].adj*chain[mu][e]
