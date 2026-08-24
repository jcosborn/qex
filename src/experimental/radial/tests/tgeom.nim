#RUNCMD $RUN1

## WP-A acceptance tests, doc/03-targets.md T1.1.

import std/[algorithm, math, sequtils, tables, unittest]
import ../core/lattice

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

const
  levs = [1, 2, 4, 8]
  tilt2 = 1.1          ## a second, deliberately different local-Lorentz gauge

func mdiff(m, n: Mat2): float =
  for i in 0..1:
    for j in 0..1: result = max(result, abs(m[i][j] - n[i][j]))

func wrap(x: float): float = x - 2.0*PI*round(x/(2.0*PI))

func dualExact(l, ls: float): float =
  ## Exact spherical replacement for the flat l l*/2: the two right spherical triangles
  ## with legs l/2 and l* have total area 4 arctan(tan(l/4) tan(l*/2)).
  4.0*arctan(tan(0.25*l)*tan(0.5*ls))

func fls(s: Sphere, f, i: int): float =
  ## l* of edge i of face f, with the sign of that face's traversal.
  s.edges[s.faces[f].e[i]].dl[if s.faces[f].s[i] > 0: 0 else: 1]

var
  sph: array[levs.len, Sphere]
  sp2: array[levs.len, Sphere]
for i, l in levs:
  sph[i] = newSphere(l)
  sp2[i] = newSphere(l, tilt2)

suite "counts and topology":
  test "N_V, N_E, N_F and Euler":                                   # T1.1a
    for i, l in levs:
      let s = sph[i]
      check s.nv == 10*l*l + 2
      check s.ne == 30*l*l
      check s.nf == 20*l*l
      check s.nv - s.ne + s.nf == 2
      check s.pos.len == s.nv
      check s.edges.len == s.ne
      check s.faces.len == s.nf

  test "every edge is in exactly two faces, with opposite orientation":
    for i, l in levs:
      let s = sph[i]
      var cnt = newSeq[int](s.ne)
      for fc in s.faces:
        for k in 0..2: inc cnt[fc.e[k]]
      check cnt.allIt(it == 2)
      check (0..<s.ne).toSeq.allIt(s.edges[it].f[0] >= 0 and s.edges[it].f[1] >= 0 and
                                   s.edges[it].f[0] != s.edges[it].f[1])
      var ok = true
      for f in 0..<s.nf:
        for k in 0..2:
          let e = s.faces[f].e[k]
          if s.faces[f].s[k] > 0:
            ok = ok and s.edges[e].f[0] == f and s.edges[e].a == s.faces[f].v[k]
          else:
            ok = ok and s.edges[e].f[1] == f and s.edges[e].b == s.faces[f].v[k]
      check ok

  test "boundary of boundary vanishes":
    for i, l in levs:
      let s = sph[i]
      var d = newSeq[int](s.nv)          # global: sum_f d(f) = 0
      var ok = true
      for f in 0..<s.nf:
        var c = initTable[int, int]()    # per face: d(d(f)) = 0
        for k in 0..2:
          let
            e = s.edges[s.faces[f].e[k]]
            h = if s.faces[f].s[k] > 0: e.b else: e.a
            t = if s.faces[f].s[k] > 0: e.a else: e.b
          c.mgetOrPut(h, 0).inc
          c.mgetOrPut(t, 0).dec
          d[h].inc
          d[t].dec
        for _, v in c: ok = ok and v == 0
      check ok
      check d.allIt(it == 0)

  test "faces are oriented outward and vertex rings are closed":
    for i, l in levs:
      let s = sph[i]
      var ok = true
      for fc in s.faces:
        let
          a = s.pos[fc.v[0]]
          b = s.pos[fc.v[1]]
          c = s.pos[fc.v[2]]
        ok = ok and dot(cross(b - a, c - a), a + b + c) > 0.0
        ok = ok and fc.area > 0.0
      check ok
      var nfc = newSeq[int](s.nv)
      for y in 0..<s.nv:
        check s.nbr[y].len == s.nbe[y].len
        check s.nbr[y].len == s.nbf[y].len
        nfc[y] = s.nbf[y].len
      var deg = newSeq[int](s.nv)
      for fc in s.faces:
        for k in 0..2: inc deg[fc.v[k]]
      check deg == nfc                   # the ring visits every incident face once

suite "areas and duals":
  test "sum A_tri = sum A_y = 4 pi":                                # T1.1b
    for i, l in levs:
      let s = sph[i]
      var sa, sy = 0.0
      for fc in s.faces: sa += fc.area
      for a in s.area: sy += a
      check abs(sa - 4.0*PI) < 1e-12
      check abs(sy - 4.0*PI) < 1e-12

  test "spherical dual decomposition of every triangle":            # T1.1c
    for i, l in levs:
      let s = sph[i]
      var r = 0.0
      for f in 0..<s.nf:
        var d = 0.0
        for k in 0..2: d += dualExact(s.edges[s.faces[f].e[k]].len, s.fls(f, k))
        r = max(r, abs(d - s.faces[f].area))
      checkpoint "L=" & $l & " max |A_tri - sum 4 atan(tan(l/4)tan(l*/2))| = " & $r
      check r < 1e-12

  test "tilde A_i partitions the triangle and the dual cell":
    for i, l in levs:
      let s = sph[i]
      var rf, ry = 0.0
      for fc in s.faces:
        rf = max(rf, abs(fc.sub[0] + fc.sub[1] + fc.sub[2] - fc.area))
      for y in 0..<s.nv:
        var t = 0.0
        for f in s.nbf[y]:
          for k in 0..2:
            if s.faces[f].v[k] == y: t += s.faces[f].sub[k]
        ry = max(ry, abs(t - s.area[y]))
      check rf < 1e-12
      check ry < 1e-12

  test "diamond area, dual length and mean edge":
    ## Edge.area is the EXACT kite area (doc/06 "THE COUPLING CONVENTION"); the
    ## flat form len*dual/2 differs at O(abar^2) and the diamonds tile the sphere.
    for i, l in levs:
      let s = sph[i]
      var r, tl, fl, ta = 0.0
      for e in s.edges:
        r = max(r, abs(e.area - dualExact(e.len, e.dl[0]) - dualExact(e.len, e.dl[1])))
        r = max(r, abs(e.dual - e.dl[0] - e.dl[1]))
        fl = max(fl, abs(e.area - 0.5*e.len*e.dual)/e.area)
        check e.area > 0.0
        check e.len > 0.0
        tl += e.len
        ta += e.area
      check r < 1e-16                    # same expression, modulo fma contraction
      check abs(ta - 4.0*PI) < 1e-12     # only the exact form tiles the sphere
      check fl < 0.1*s.abar*s.abar       # flat form is O(abar^2) away
      check fl > 1e-4*s.abar*s.abar
      check abs(tl - float(s.ne)*s.abar) < 1e-12

  test "abar table":                                                # T1.1h
    check abs(sph[0].abar - 2.0*arcsin(sqrt(0.1*(5.0 - sqrt 5.0)))) < 1e-14
    for i, l in levs:
      checkpoint "L=" & $l & " abar = " & $sph[i].abar
    check abs(sph[0].abar - 1.107148718) < 5e-10
    check abs(sph[1].abar - 0.590946445) < 5e-10
    check abs(sph[2].abar - 0.299474473) < 5e-10
    check abs(sph[3].abar - 0.150227491) < 5e-10

  test "flat form of the dual decomposition is only O(abar^2)":
    ## sum_i l_i l*_i/2 = A_tri holds exactly in the plane, not on the sphere.
    ## The relative residual must fall like abar^2.
    var rel: array[levs.len, float]
    for i, l in levs:
      let s = sph[i]
      var r = 0.0
      for f in 0..<s.nf:
        var d = 0.0
        for k in 0..2: d += 0.5*s.edges[s.faces[f].e[k]].len*s.fls(f, k)
        r = max(r, abs(d - s.faces[f].area))
      rel[i] = r*float(s.nf)/(4.0*PI)
      checkpoint "L=" & $l & " relative residual = " & $rel[i]
    for i in 1..<levs.len: check rel[i] < rel[i-1]
    for i in 2..<levs.len: check rel[i-1]/rel[i] > 2.5

suite "simplicial closure relation (IV.6)":                          # T1.1d
  test "off-diagonal and diagonal residuals are O(abar^2)":
    var od, dg: array[levs.len, float]
    for i, l in levs:
      let
        s = sph[i]
        c = s.checkClosure
        a3 = 4.0*PI/float(s.nf)
      od[i] = c.offDiag/a3
      dg[i] = c.diag/a3
      checkpoint "L=" & $l & " offDiag/A_tri = " & $od[i] & "  |diag-A_tri|/A_tri = " & $dg[i]
    check od[0] < 1e-14                  # L=1 is exact by the icosahedral 3-fold symmetry
    for i in 1..<levs.len: check dg[i] < dg[i-1]
    for i in 2..<levs.len:
      check dg[i-1]/dg[i] > 2.5
      check od[i-1]/od[i] > 2.5

suite "spin connection":
  test "face holonomy = exp(i sigma3 A_tri/2)":                      # T1.1e
    for i, l in levs:
      var r = 0.0
      for s in [sph[i], sp2[i]]:
        for f in 0..<s.nf: r = max(r, mdiff(s.holonomy f, expIsig3(s.faces[f].area)))
      checkpoint "L=" & $l & " max |holonomy - expIsig3(A_tri)| = " & $r
      check r < 1e-12

  test "global holonomy is the identity":                            # T1.1f
    for i, l in levs:
      let s = sph[i]
      var g = id2
      for f in 0..<s.nf: g = g*s.holonomy(f)
      checkpoint "L=" & $l & " |prod holonomy - 1| = " & $mdiff(g, id2)
      check mdiff(g, id2) < 1e-11

  test "primary omega agrees with the omegaChart oracle mod 2 pi":
    for i, l in levs:
      let s = sph[i]
      var r, rq = 0.0
      for e in 0..<s.ne:
        r = max(r, abs(wrap(s.omegaChart(e) - s.edges[e].omega)))
        rq = max(rq, abs(s.omegaChart(e, 64) - s.omegaChart(e, 32)))
      checkpoint "L=" & $l & " max |omega - omegaChart| mod 2pi = " & $r
      check r < 1e-12
      check rq < 1e-12                   # the oracle unwraps exactly, so nq must not matter

  test "the F2 lift minus the chart lift is exactly the pole cut":
    ## The raw chart lift expIsig3(omegaChart) is not a spin structure: by (B.4) its
    ## holonomy is -exp(i sigma3 A/2) on the two faces that contain the chart poles.
    ## The ratio r = sgn*(-1)^winding must therefore have coboundary -1 on exactly those
    ## two faces.  Cancelling that with a dual path between them -- the discrete meridian
    ## cut of doc/02 section 2.3 -- leaves a pure vertex coboundary.
    for i, l in levs:
      let s = sph[i]
      var r = newSeq[float](s.ne)
      for e in 0..<s.ne:
        let k = int round((s.omegaChart(e) - s.edges[e].omega)/(2.0*PI))
        r[e] = s.edges[e].sgn*(if (k and 1) == 1: -1.0 else: 1.0)
        check abs(abs(s.edges[e].sgn) - 1.0) < 1e-30
      var pole = [-1, -1]
      for f in 0..<s.nf:
        let fc = s.faces[f]
        for j in 0..1:
          let z = if j == 0: s.chart[2] else: -1.0*s.chart[2]
          if sphArea(s.pos[fc.v[0]], s.pos[fc.v[1]], z) > 0.0 and
             sphArea(s.pos[fc.v[1]], s.pos[fc.v[2]], z) > 0.0 and
             sphArea(s.pos[fc.v[2]], s.pos[fc.v[0]], z) > 0.0: pole[j] = f
      check pole[0] >= 0
      check pole[1] >= 0
      check pole[0] != pole[1]

      # sum_f sum_i s_i omegaChart_i = 0 (each edge twice, opposite signs) while
      # sum_f A_tri = 4 pi, so the chart lift is short by exactly 2 pi on each of
      # chi(S^2) = 2 faces -- and 2 pi is a sign, not nothing: expIsig3 has period 4 pi.
      var
        dtot = 0.0
        ndef = 0
      for f in 0..<s.nf:
        var d = -s.faces[f].area
        for k in 0..2: d += float(s.faces[f].s[k])*s.omegaChart(s.faces[f].e[k])
        dtot += d
        if abs(d) > PI:
          inc ndef
          check (f == pole[0] or f == pole[1])
          check abs(d + 2.0*PI) < 1e-12
      check ndef == 2
      check abs(dtot + 4.0*PI) < 1e-11
      var nbad = 0
      for f in 0..<s.nf:
        if r[s.faces[f].e[0]]*r[s.faces[f].e[1]]*r[s.faces[f].e[2]] < 0.0:
          inc nbad
          check (f == pole[0] or f == pole[1])
      check nbad == 2

      # the cut: a dual path pole[0] -> pole[1]; every edge it crosses flips
      var
        par = newSeq[int](s.nf)
        pare = newSeq[int](s.nf)
        qf = @[pole[0]]
        hf = 0
      for f in 0..<s.nf: par[f] = -1
      par[pole[0]] = pole[0]
      while hf < qf.len:
        let f = qf[hf]
        inc hf
        for k in 0..2:
          let
            e = s.faces[f].e[k]
            g = if s.edges[e].f[0] == f: s.edges[e].f[1] else: s.edges[e].f[0]
          if par[g] < 0:
            par[g] = f
            pare[g] = e
            qf.add g
      var fc = pole[1]
      while fc != pole[0]:
        r[pare[fc]] = -r[pare[fc]]
        fc = par[fc]
      for f in 0..<s.nf:
        check r[s.faces[f].e[0]]*r[s.faces[f].e[1]]*r[s.faces[f].e[2]] > 0.0

      var
        tau = newSeq[float](s.nv)
        seen = newSeq[bool](s.nv)
        q = @[0]
        h = 0
      tau[0] = 1.0
      seen[0] = true
      while h < q.len:
        let y = q[h]
        inc h
        for k, n in s.nbr[y]:
          if not seen[n]:
            seen[n] = true
            tau[n] = r[s.nbe[y][k]]*tau[y]
            q.add n
      check seen.allIt(it)
      var bad = 0
      for e in 0..<s.ne:
        if r[e]*tau[s.edges[e].a]*tau[s.edges[e].b] < 0.0: inc bad
      check bad == 0

  test "eb is ea rotated by omega":
    for i, l in levs:
      let s = sph[i]
      var r = 0.0
      for e in s.edges:
        let
          c = cos e.omega
          w = sin e.omega
        r = max(r, abs(e.eb[0] - (c*e.ea[0] - w*e.ea[1])))
        r = max(r, abs(e.eb[1] - (w*e.ea[0] + c*e.ea[1])))
        r = max(r, abs(e.ea[0]*e.ea[0] + e.ea[1]*e.ea[1] - 1.0))
      check r < 1e-13

  test "parity: omega(P y1, P y2) = -omega(y1, y2)":                 # (III.3), (IV.16)
    ## P: (theta, phi) -> (pi-theta, phi+pi) is the antipodal map in any polar chart,
    ## and the icosahedron is centrally symmetric, so this is an exact lattice identity.
    for i, l in levs:
      let s = sph[i]
      var p = newSeq[int](s.nv)
      for y in 0..<s.nv:
        p[y] = -1
        for z in 0..<s.nv:
          let d = s.pos[y] + s.pos[z]
          if dot(d, d) < 1e-18: p[y] = z
      check p.allIt(it >= 0)
      var emap = initTable[(int, int), int]()
      for e in 0..<s.ne: emap[(s.edges[e].a, s.edges[e].b)] = e
      var r, ra = 0.0
      for e in 0..<s.ne:
        let
          pa = p[s.edges[e].a]
          pb = p[s.edges[e].b]
          fwd = emap.getOrDefault((min(pa, pb), max(pa, pb)), -1)
        check fwd >= 0
        # canonical a<b may swap the traversal, which flips omega a second time
        let sg = if pa < pb: -1.0 else: 1.0
        r = max(r, abs(s.edges[fwd].omega - sg*s.edges[e].omega))
        ra = max(ra, abs(s.edges[fwd].len - s.edges[e].len))
      for y in 0..<s.nv: ra = max(ra, abs(s.area[p[y]] - s.area[y]))
      checkpoint "L=" & $l & " max |omega(Pe) + omega(e)| = " & $r
      check r < 1e-13
      check ra < 1e-14

  test "chart poles are clear of every site and link geodesic":
    for i, l in levs:
      let
        s = sph[i]
        g = s.poleGap
      checkpoint "L=" & $l & " poleGap site/link = " & $g.site & " / " & $g.link
      check g.site > 0.2*s.abar
      check g.link > 0.1*s.abar

suite "tilt independence":
  test "metric data does not involve the chart at all":
    for i, l in levs:
      let
        s = sph[i]
        t = sp2[i]
      check s.abar == t.abar
      check s.area == t.area
      for e in 0..<s.ne:
        check s.edges[e].len == t.edges[e].len
        check s.edges[e].dl == t.edges[e].dl
        check s.edges[e].area == t.edges[e].area
      for f in 0..<s.nf:
        check s.faces[f].area == t.faces[f].area
        check s.faces[f].sub == t.faces[f].sub

  test "face holonomies are gauge invariant":
    for i, l in levs:
      let
        s = sph[i]
        t = sp2[i]
      var r, d = 0.0
      for f in 0..<s.nf: r = max(r, mdiff(s.holonomy f, t.holonomy f))
      for e in 0..<s.ne: d = max(d, abs(s.edges[e].omega - t.edges[e].omega))
      check r < 1e-12
      check d > 0.1                      # the gauge really is different

suite "flat equilateral cross-check":
  ## Regular hexagonal patch of polar radius d about the north pole.  The six triangles
  ## become equilateral as d -> 0, so l* -> a/(2 sqrt3), A_tri -> sqrt3 a^2/4,
  ## A_y -> sqrt3 a^2/2, kappa -> 1/sqrt3, kappa' -> (sqrt3/2)(a/at), all to O(d^2).
  proc patch(d: float): tuple[ls, atri, ay, kap, kapT: float] =
    let o: Vec3 = [0.0, 0.0, 1.0]
    var p, cc: array[6, Vec3]
    for k in 0..5:
      let ph = float(k)*PI/3.0
      p[k] = [sin(d)*cos(ph), sin(d)*sin(ph), cos(d)]
    for k in 0..5: cc[k] = circum(o, p[k], p[(k+1) mod 6])
    var a = 0.0
    for k in 0..5: a += geodesic(o, p[k]) + geodesic(p[k], p[(k+1) mod 6])
    a /= 12.0
    let
      l0 = dualLen(cc[0], o, p[0])
      l1 = dualLen(cc[5], p[0], o)
    result.ls = 0.5*(l0 + l1)/a
    result.atri = sphArea(o, p[0], p[1])/(a*a)
    for k in 0..5: result.ay += sphArea(o, cc[k], cc[(k+1) mod 6])
    result.ay /= a*a
    result.kap = (l0 + l1)/a
    result.kapT = result.ay                 # kappa'/(a/at) = A_y/a^2

  test "flat relations hold with an O(d^2) error":
    let
      r3 = sqrt 3.0
      c = patch 1e-2
      f = patch 1e-3
      ec = [abs(c.ls - 0.5/r3), abs(c.atri - 0.25*r3), abs(c.ay - 0.5*r3),
            abs(c.kap - 1.0/r3), abs(c.kapT - 0.5*r3)]
      ef = [abs(f.ls - 0.5/r3), abs(f.atri - 0.25*r3), abs(f.ay - 0.5*r3),
            abs(f.kap - 1.0/r3), abs(f.kapT - 0.5*r3)]
    checkpoint "d=1e-2 errors " & $ec
    checkpoint "d=1e-3 errors " & $ef
    for k in 0..4:
      check ec[k] < 2e-5                    # O(d^2) at d = 1e-2
      check ef[k] < 2e-7                    # 100x smaller at d = 1e-3
      check ef[k] > 1e-12                   # and it really is the d^2 term, not roundoff

  test "mean kappa on the sphere tends to 1/sqrt3":
    var e: array[levs.len, float]
    for i, l in levs:
      let lat = newLat(sph[i], 8, 0.2)
      var k = 0.0
      for x in lat.kap: k += x
      e[i] = abs(k/float(sph[i].ne) - 1.0/sqrt(3.0))
      checkpoint "L=" & $l & " |<kappa> - 1/sqrt3| = " & $e[i]
    for i in 1..<levs.len: check e[i-1]/e[i] > 3.5

suite "icosahedral orbits":                                          # T1.1g
  test "12 five-fold sites, the rest six-fold":
    for i, l in levs:
      let s = sph[i]
      var n5, n6 = 0
      for y in 0..<s.nv:
        if s.nbr[y].len == 5: inc n5
        elif s.nbr[y].len == 6: inc n6
      check n5 == 12
      check n6 == s.nv - 12

  test "A_y takes one value per I_h orbit":
    ## The orbits are the barycentric multisets {i,j,k}, i+j+k = L, so their number is
    ## the count of partitions of L into at most three parts, round((L+3)^2/12).
    for i, l in levs:
      let s = sph[i]
      var a = s.area
      a.sort
      var n = 1
      for k in 1..<a.len:
        if a[k] - a[k-1] > 1e-12: inc n
      checkpoint "L=" & $l & " distinct A_y = " & $n
      check n == int round(float((l+3)*(l+3))/12.0)

suite "lattice couplings":
  test "kappa, kappa', volume weights":
    let
      s = sph[2]
      at = 0.2
      lat = newLat(s, 12, at)
    check lat.nsite == s.nv*12
    var r = 0.0
    for e in 0..<s.ne:
      # the flat identity kap = dual/abar holds only to O(abar^2) now
      r = max(r, abs(lat.kap[e] - 2.0*s.edges[e].area/(s.abar*s.edges[e].len)))
    var v = 0.0
    for y in 0..<s.nv:
      r = max(r, abs(lat.kapT[y] - s.area[y]/(s.abar*at)))
      r = max(r, abs(lat.volw[y] - s.area[y]*at))
      v += lat.volw[y]
    check r < 1e-15
    check abs(lat.volbar - v/float(s.nv)) < 1e-15
    check abs(lat.volbar - 4.0*PI*at/float(s.nv)) < 1e-14

  test "gauge couplings and the doubler window":
    let
      s = sph[2]
      lat = newLat(s, 12, 0.2)
      g2 = 1.7
    check abs(lat.betaFace(5, g2) - lat.at/(g2*s.faces[5].area)) < 1e-15
    check abs(lat.betaEdge(7, g2) -
              2.0*s.edges[7].area/(g2*s.edges[7].len*s.edges[7].len*lat.at)) < 1e-15
    check abs(lat.asOverAt - s.abar/lat.at) < 1e-15
    check lat.asOverAt >= 4.0/3.0                       # L=4, at=0.2 still passes
    check abs(lat.maxM - 0.9*min(4.0/sqrt(3.0), sqrt(3.0)*lat.asOverAt)) < 1e-15
    check newLat(sph[3], 12, 0.2).asOverAt < 4.0/3.0    # L=8 at at=0.2 does not

  test "kappa and A_y ranges match the independent oracle":
    ## Pinned against the pure-Python geometry oracle of doc/06-status.md.
    const
      # exact-kappa oracle values (doc/06 "THE COUPLING CONVENTION" section);
      # the flat-kappa values were 0.65911 / 0.53248-0.66331 / 0.45471-0.61252 / 0.37082-0.62676
      kr = [[0.683450, 0.683450], [0.538115, 0.669683], [0.455754, 0.613770],
            [0.370994, 0.627216]]
      ar = [[1.047198, 1.047198], [0.273844, 0.309341], [0.058328, 0.083977],
            [0.012966, 0.022939]]
    for i, l in levs:
      let
        s = sph[i]
        lat = newLat(s, 8, 0.2)
      var
        kmin = lat.kap[0]
        kmax = lat.kap[0]
        amin = s.area[0]
        amax = s.area[0]
      for k in lat.kap:
        kmin = min(kmin, k)
        kmax = max(kmax, k)
      for a in s.area:
        amin = min(amin, a)
        amax = max(amax, a)
      check abs(kmin - kr[i][0]) < 5e-6
      check abs(kmax - kr[i][1]) < 5e-6
      check abs(amin - ar[i][0]) < 5e-7
      check abs(amax - ar[i][1]) < 5e-7

  test "Table I doubler membership, T=16":                           # T1.4i
    ## abar/at >= 4/3 with at = T/Lt is Lt >= T/(0.75 abar); at T = 16 this reproduces
    ## the paper's Table I, all 20 cells.
    const
      lts = [64, 96, 120, 144, 168]
      want = [[true, true, true, true, true],      # L=1, Lt_min = 19.3
              [true, true, true, true, true],      # L=2, Lt_min = 36.1
              [false, true, true, true, true],     # L=4, Lt_min = 71.2
              [false, false, false, true, true]]   # L=8, Lt_min = 142.0
      lmin = [19.3, 36.1, 71.2, 142.0]
      tex = 16.0
    for i, l in levs:
      let s = sph[i]
      checkpoint "L=" & $l & " Lt_min = " & $(tex/(0.75*s.abar))
      check abs(tex/(0.75*s.abar) - lmin[i]) < 0.1
      for j, lt in lts:
        let lat = newLat(s, lt, tex/float(lt))
        check (lat.asOverAt >= 4.0/3.0) == want[i][j]

  test "anisotropy at the interacting-run a_t = 0.2":
    const ao = [5.54, 2.95, 1.50, 0.75]
    for i, l in levs:
      let lat = newLat(sph[i], 8, 0.2)
      checkpoint "L=" & $l & " abar/at = " & $lat.asOverAt
      check abs(lat.asOverAt - ao[i]) < 5e-3
      check (lat.asOverAt >= 4.0/3.0) == (l <= 4)

  test "index maps wrap in t":
    let lat = newLat(sph[0], 6, 0.2)
    check lat.sIdx(3, 0) == 3
    check lat.sIdx(3, 6) == 3
    check lat.sIdx(3, -1) == lat.sIdx(3, 5)
    check lat.eIdx(2, -1) == 2 + sph[0].ne*5
    check lat.tIdx(2, 7) == 2 + sph[0].nv
