## Refined icosahedron on S^2: sites, links, triangles, their circumcentric duals,
## and the lattice spin connection.  Radius R = 1, so every length is an arc length
## and every area a spherical area -- the complex is intrinsic, not the flat polyhedron.
##
## Normative reference: doc/02-formulation.md section 2, doc/04-interfaces.md section 3.

import std/[math, tables]
import types
export types

const
  defaultTilt* = 0.7
    ## Polar tilt of the local-Lorentz chart.  A gauge choice: it must keep every site
    ## and every link geodesic clear of the chart poles (see `poleGap`), and no
    ## observable may depend on it.  0.7 was picked by scanning `poleGap` over L = 1..16:
    ## it leaves the poles ~abar/2 from the nearest site and ~abar/4 from the nearest
    ## link geodesic at every level, which is about the best a single axis can do.

type
  Edge* = object
    a*, b*: int              ## endpoints, canonical orientation a -> b, a < b
    len*: float              ## geodesic length l_e
    dl*: array[2, float]     ## signed dual lengths l*, one per incident face, ordered as `f`
    dual*: float             ## dl[0] + dl[1]
    area*: float             ## diamond area A_e, EXACT spherical (kiteArea sum over dl);
                             ## the flat form is len*dual/2, O(abar^2) apart
    ea*: array[2, float]     ## e^a_{ab}(a) = (cos alpha_a, sin alpha_a)
    eb*: array[2, float]     ## e^a_{ab}(b): the tangent at b continuing away from a,
                             ## so the reverse hop needs e^a_{ba}(b) = -eb
    omega*: float            ## spin-connection angle, principal value; Omega = sgn*expIsig3(omega)
    sgn*: float              ## +-1, the Z2 spin-structure sign
    f*: array[2, int]        ## incident faces; f[0] traverses a->b, f[1] traverses b->a

  Face* = object
    v*: array[3, int]        ## vertices, counterclockwise seen from outside
    e*: array[3, int]        ## e[i] joins v[i] and v[(i+1) mod 3]
    s*: array[3, int]        ## +1 if e[i] runs v[i]->v[i+1], else -1
    area*: float             ## spherical area A_tri (the excess)
    cc*: Vec3                ## dual point: circumcenter on S^2
    sub*: array[3, float]    ## tilde A_i for v[0..2]; sum = area

  Sphere* = ref object
    lev*: int                ## refinement level L
    nv*, ne*, nf*: int       ## 10L^2+2, 30L^2, 20L^2
    pos*: seq[Vec3]          ## unit position vectors
    area*: seq[float]        ## A_y, dual polygon area per site
    nbr*: seq[seq[int]]      ## neighbour sites, counterclockwise around the site
    nbe*: seq[seq[int]]      ## nbe[y][k] is the edge joining y and nbr[y][k]
    nbf*: seq[seq[int]]      ## nbf[y][k] is the face (y, nbr[y][k], nbr[y][k+1])
    edges*: seq[Edge]
    faces*: seq[Face]
    abar*: float             ## mean edge length
    tilt*: float
    chart*: array[3, Vec3]   ## chart axes (x', y', z'); z' is the polar axis

# --- spherical primitives ---------------------------------------------------

func sphArea*(a, b, c: Vec3): float =
  ## Signed spherical excess (van Oosterom-Strackee):
  ##   tan(A/2) = a.((b-a)x(c-a)) / (1 + a.b + b.c + c.a).
  ## Positive for a,b,c counterclockwise seen from outside.  The (b-a)x(c-a) form
  ## keeps the numerator free of the cancellation of a.(bxc) for small triangles.
  2.0*arctan2(dot(a, cross(b - a, c - a)), 1.0 + dot(a, b) + dot(b, c) + dot(c, a))

func circum*(a, b, c: Vec3): Vec3 =
  ## Dual point: circumcenter of the spherical triangle, on the a+b+c hemisphere.
  let n = unit cross(b - a, c - a)
  if dot(n, a + b + c) > 0.0: n else: -1.0*n

func dualLen*(cc, p, q: Vec3): float =
  ## Signed dual length l* = arcsin(cc . unit(p x q)) for the edge p->q,
  ## positive when cc lies on the interior side of the edge.
  arcsin dot(cc, unit cross(p, q))

func kiteArea*(le, ls: float): float =
  ## Exact spherical area of the kite of an edge on one side: the two right
  ## triangles with legs le/2 and ls,  E = 4 arctan(tan(le/4) tan(ls/2)).
  ## Odd in ls, so signed dual lengths work.  Sums to A_tri over a face exactly;
  ## the flat form le*ls/2 does so only to O(abar^2).  This is the paper's A_l:
  ## both published Delta_0 values pin it (doc/06 "THE COUPLING CONVENTION").
  4.0*arctan(tan(0.25*le)*tan(0.5*ls))

func ptrans*(p, q, t: Vec3): Vec3 =
  ## Parallel transport of the tangent t from p to q along the geodesic.
  ## gamma(s) = cos(s) p + sin(s) d_p, so d_q = unit(cos(th) q - p); the component
  ## along n = unit(p x q) is untouched.
  let
    n = unit cross(p, q)
    c = dot(p, q)
    dp = unit(q - c*p)
    dq = unit(c*q - p)
  dot(t, dp)*dq + dot(t, n)*n

func chartFrame*(tilt: float): array[3, Vec3] =
  ## Axes of the polar chart.  z' is tilted by `tilt` off the icosahedral z axis at
  ## azimuth phi*tilt, an irrational multiple that misses every icosahedral symmetry axis.
  let
    ct = cos tilt
    st = sin tilt
    ph = 1.6180339887498949*tilt
    cp = cos ph
    sp = sin ph
    z = [st*cp, st*sp, ct]
    x = [ct*cp, ct*sp, -st]
  [x, cross(z, x), z]

func pframe*(ch: array[3, Vec3], p: Vec3): array[2, Vec3] =
  ## (e_theta, e_phi) of the chart `ch` at p.  e_theta x e_phi = p, so the frame is
  ## oriented by the outward normal everywhere.
  let
    ct = dot(ch[2], p)
    x = dot(ch[0], p)
    y = dot(ch[1], p)
    st = sqrt(x*x + y*y)
    cp = x/st
    sp = y/st
  [ct*(cp*ch[0] + sp*ch[1]) - st*ch[2], (-sp)*ch[0] + cp*ch[1]]

# --- construction -----------------------------------------------------------

func icosa(): tuple[v: seq[Vec3], f: seq[array[3, int]]] =
  ## The 12 cyclic permutations of (0, +-1, +-phi)/sqrt(1+phi^2) and the 20 faces,
  ## each oriented outward.
  let
    p = 0.5*(1.0 + sqrt 5.0)
    s = 1.0/sqrt(1.0 + p*p)
  for i in 0..3:
    let
      u = if (i and 1) == 0: s else: -s
      w = if (i and 2) == 0: p*s else: -p*s
    result.v.add [0.0, u, w]
    result.v.add [u, w, 0.0]
    result.v.add [w, 0.0, u]
  let cut = 6.0*s*s          # edge^2 = 4 s^2, next chord^2 = 4 p^2 s^2 = 10.5 s^2
  for i in 0..11:
    for j in i+1..11:
      let dij = result.v[i] - result.v[j]
      if dot(dij, dij) > cut: continue
      for k in j+1..11:
        let
          dik = result.v[i] - result.v[k]
          djk = result.v[j] - result.v[k]
        if dot(dik, dik) > cut or dot(djk, djk) > cut: continue
        let (a, b, c) = (result.v[i], result.v[j], result.v[k])
        if dot(cross(b - a, c - a), a + b + c) > 0.0: result.f.add [i, j, k]
        else: result.f.add [i, k, j]

proc solveZ2(nrow, ncol: int, rows: var seq[seq[uint64]], rhs: var seq[bool]): seq[bool] =
  ## Gaussian elimination over F2, free variables set to 0.  For the face-edge
  ## incidence matrix of S^2 the rank is nrow-1 and the system is consistent
  ## because sum_f k_f = 0.
  let nw = (ncol + 63) div 64
  var
    piv: seq[int]
    r = 0
  for c in 0..<ncol:
    let
      w = c shr 6
      m = 1'u64 shl (c and 63)
    var p = -1
    for i in r..<nrow:
      if (rows[i][w] and m) != 0:
        p = i
        break
    if p < 0: continue
    swap(rows[r], rows[p])
    swap(rhs[r], rhs[p])
    for i in 0..<nrow:
      if i != r and (rows[i][w] and m) != 0:
        for k in 0..<nw: rows[i][k] = rows[i][k] xor rows[r][k]
        rhs[i] = rhs[i] xor rhs[r]
    piv.add c
    inc r
  result = newSeq[bool](ncol)
  for i in 0..<r: result[piv[i]] = rhs[i]

proc newSphere*(lev: int, tilt = defaultTilt): Sphere =
  ## Level-`lev` refined icosahedron projected onto S^2.  `tilt` selects the
  ## local-Lorentz gauge (doc/02 section 2.3); observables must not depend on it.
  let ico = icosa()
  var
    pos: seq[Vec3]
    grid = initTable[array[3, int64], seq[int]]()
  proc vidx(p: Vec3): int =
    ## Deduplicate shared refinement vertices; coordinates are hashed on a 1e-9 grid
    ## and the 26 neighbouring cells are probed so a point never straddles a boundary.
    let k = [int64 round(p[0]*1e9), int64 round(p[1]*1e9), int64 round(p[2]*1e9)]
    for dx in -1'i64..1'i64:
      for dy in -1'i64..1'i64:
        for dz in -1'i64..1'i64:
          let kk = [k[0]+dx, k[1]+dy, k[2]+dz]
          if grid.hasKey kk:
            for i in grid[kk]:
              let d = pos[i] - p
              if dot(d, d) < 1e-18: return i
    result = pos.len
    pos.add p
    grid.mgetOrPut(k, @[]).add result

  # v_ijk = unit(i A + j B + k C), i+j+k = L; L^2 subtriangles per icosahedral face
  var tri: seq[array[3, int]]
  for fc in ico.f:
    let
      a = ico.v[fc[0]]
      b = ico.v[fc[1]]
      c = ico.v[fc[2]]
    var id = newSeq[seq[int]](lev+1)
    for i in 0..lev:
      id[i] = newSeq[int](lev+1-i)
      for j in 0..lev-i:
        id[i][j] = vidx unit(float(i)*a + float(j)*b + float(lev-i-j)*c)
    for i in 0..<lev:
      for j in 0..<lev-i:
        tri.add [id[i+1][j], id[i][j+1], id[i][j]]
        if i+j < lev-1: tri.add [id[i+1][j], id[i][j+1], id[i+1][j+1]]

  result = Sphere(lev: lev, tilt: tilt, chart: chartFrame(tilt), pos: pos,
                  nv: pos.len, nf: tri.len)
  result.faces = newSeq[Face](result.nf)
  var emap = initTable[(int, int), int]()
  for n, t0 in tri:
    var t = t0
    if dot(cross(pos[t[1]] - pos[t[0]], pos[t[2]] - pos[t[0]]),
           pos[t[0]] + pos[t[1]] + pos[t[2]]) < 0.0: swap(t[1], t[2])
    result.faces[n].v = t
    for i in 0..2:
      let
        p = t[i]
        q = t[(i+1) mod 3]
        k = (min(p, q), max(p, q))
      var e = emap.getOrDefault(k, -1)
      if e < 0:
        e = result.edges.len
        emap[k] = e
        result.edges.add Edge(a: k[0], b: k[1], f: [-1, -1])
      result.faces[n].e[i] = e
      if p == k[0]:
        result.faces[n].s[i] = 1
        result.edges[e].f[0] = n
      else:
        result.faces[n].s[i] = -1
        result.edges[e].f[1] = n
  result.ne = result.edges.len

  for n in 0..<result.nf:
    let
      v = result.faces[n].v
      a = pos[v[0]]
      b = pos[v[1]]
      c = pos[v[2]]
      cc = circum(a, b, c)
      m = [unit(a + b), unit(b + c), unit(c + a)]
    result.faces[n].cc = cc
    result.faces[n].area = sphArea(a, b, c)
    # tilde A_i = quadrilateral (v_i, m_i, cc, m_{i-1}), split at cc
    for i in 0..2:
      let p = pos[v[i]]
      result.faces[n].sub[i] = sphArea(p, m[i], cc) + sphArea(p, cc, m[(i+2) mod 3])

  var tot = 0.0
  for e in 0..<result.ne:
    let
      a = pos[result.edges[e].a]
      b = pos[result.edges[e].b]
    result.edges[e].len = geodesic(a, b)
    result.edges[e].dl[0] = dualLen(result.faces[result.edges[e].f[0]].cc, a, b)
    result.edges[e].dl[1] = dualLen(result.faces[result.edges[e].f[1]].cc, b, a)
    result.edges[e].dual = result.edges[e].dl[0] + result.edges[e].dl[1]
    result.edges[e].area = kiteArea(result.edges[e].len, result.edges[e].dl[0]) +
                           kiteArea(result.edges[e].len, result.edges[e].dl[1])
    tot += result.edges[e].len
  result.abar = tot/float(result.ne)

  # counterclockwise vertex rings: in a face (v0,v1,v2) the neighbour after v1 is v2
  var nxt = initTable[(int, int), (int, int)]()
  for n in 0..<result.nf:
    let v = result.faces[n].v
    for i in 0..2: nxt[(v[i], v[(i+1) mod 3])] = (v[(i+2) mod 3], n)
  var first = newSeq[int](result.nv)
  for y in 0..<result.nv: first[y] = result.nv
  for e in 0..<result.ne:
    let
      a = result.edges[e].a
      b = result.edges[e].b
    first[a] = min(first[a], b)
    first[b] = min(first[b], a)
  result.nbr = newSeq[seq[int]](result.nv)
  result.nbe = newSeq[seq[int]](result.nv)
  result.nbf = newSeq[seq[int]](result.nv)
  result.area = newSeq[float](result.nv)
  for y in 0..<result.nv:
    var p = first[y]
    while true:
      let (q, fc) = nxt[(y, p)]
      result.nbr[y].add p
      result.nbe[y].add emap[(min(y, p), max(y, p))]
      result.nbf[y].add fc
      p = q
      if p == first[y]: break
    let n = result.nbf[y].len
    for k in 0..<n:
      result.area[y] += sphArea(pos[y], result.faces[result.nbf[y][k]].cc,
                                result.faces[result.nbf[y][(k+1) mod n]].cc)

  # tangent-frame components and the primary spin connection: transport e_theta(a)
  # along the geodesic and read its angle in the frame at b.
  for e in 0..<result.ne:
    let
      a = pos[result.edges[e].a]
      b = pos[result.edges[e].b]
      fa = pframe(result.chart, a)
      fb = pframe(result.chart, b)
      ta = tangent(a, b)
      tb = ptrans(a, b, ta)
      pt = ptrans(a, b, fa[0])
    result.edges[e].ea = [dot(ta, fa[0]), dot(ta, fa[1])]
    result.edges[e].eb = [dot(tb, fb[0]), dot(tb, fb[1])]
    result.edges[e].omega = arctan2(dot(pt, fb[1]), dot(pt, fb[0]))

  # Z2 lift: prod_i sgn_i = (-1)^k_f with sum_i s_i omega_i - A_tri = 2 pi k_f.
  let nw = (result.ne + 63) div 64
  var
    rows = newSeq[seq[uint64]](result.nf)
    rhs = newSeq[bool](result.nf)
  for n in 0..<result.nf:
    rows[n] = newSeq[uint64](nw)
    var w = -result.faces[n].area
    for i in 0..2:
      let e = result.faces[n].e[i]
      w += float(result.faces[n].s[i])*result.edges[e].omega
      rows[n][e shr 6] = rows[n][e shr 6] xor (1'u64 shl (e and 63))
    rhs[n] = (int(round(w/(2.0*PI))) and 1) == 1
  let x = solveZ2(result.nf, result.ne, rows, rhs)
  for e in 0..<result.ne:
    result.edges[e].sgn = if x[e]: -1.0 else: 1.0

# --- oracle, holonomy, diagnostics ------------------------------------------

proc omegaChart*(s: Sphere, e: int, nq = 32): float =
  ## Independent oracle for (III.2), omega = -int cos(theta) dphi.  The chart angle
  ## alpha(s) of the geodesic tangent is sampled and unwrapped, so the returned value
  ## is alpha(l) - alpha(0) exactly -- not a quadrature -- provided no step wraps by pi.
  let
    p = s.pos[s.edges[e].a]
    q = s.pos[s.edges[e].b]
    th = geodesic(p, q)
    d = tangent(p, q)
  var prev = 0.0
  for k in 0..nq:
    let
      u = th*float(k)/float(nq)
      x = cos(u)*p + sin(u)*d
      t = (-sin(u))*p + cos(u)*d
      fr = pframe(s.chart, x)
      al = arctan2(dot(t, fr[1]), dot(t, fr[0]))
    if k > 0:
      var dd = al - prev
      while dd > PI: dd -= 2.0*PI
      while dd < -PI: dd += 2.0*PI
      result += dd
    prev = al

proc holonomy*(s: Sphere, f: int): Mat2 =
  ## Product of Omega around face f in the orientation of Face.s.  Must equal
  ## expIsig3(A_tri) -- doc/02 (B.4), slide 7.
  result = id2
  for i in 0..2:
    let e = s.faces[f].e[i]
    result = result*(cr(s.edges[e].sgn)*expIsig3(float(s.faces[f].s[i])*s.edges[e].omega))

proc checkClosure*(s: Sphere): tuple[offDiag, diag: float] =
  ## Simplicial closure relation (IV.6): sum_i l_i l*_i e^a_i e^b_i = A_tri delta^ab
  ## up to O(abar^2).  The edge tangents are taken at the edge midpoints and parallel
  ## transported to the dual point, where all three are read in one frame.
  for f in 0..<s.nf:
    let
      fc = s.faces[f]
      c = fc.cc
      p0 = s.pos[fc.v[0]]
      u = unit(p0 - dot(c, p0)*c)
      w = cross(c, u)
    var m: array[2, array[2, float]]
    for i in 0..2:
      let
        p = s.pos[fc.v[i]]
        q = s.pos[fc.v[(i+1) mod 3]]
        mid = unit(p + q)
        t = ptrans(mid, c, tangent(mid, q))
        e0 = dot(t, u)
        e1 = dot(t, w)
        g = s.edges[fc.e[i]].len*s.edges[fc.e[i]].dl[if fc.s[i] > 0: 0 else: 1]
      m[0][0] += g*e0*e0
      m[0][1] += g*e0*e1
      m[1][1] += g*e1*e1
    result.offDiag = max(result.offDiag, abs(m[0][1]))
    result.diag = max(result.diag, max(abs(m[0][0] - fc.area), abs(m[1][1] - fc.area)))

proc poleGap*(s: Sphere): tuple[site, link: float] =
  ## Angular distance from the chart poles to the nearest site and to the nearest point
  ## of any link geodesic.  Both must stay well away from 0: sites need a frame, and
  ## `omegaChart` needs alpha to move by less than pi between samples.
  ## On the segment x(u) = cos(u) p + sin(u) d one has x.z = R cos(u-u0), R = hypot(p.z, d.z).
  result.site = PI
  result.link = PI
  for p in s.pos:
    result.site = min(result.site, arccos(min(1.0, abs(dot(s.chart[2], p)))))
  for e in s.edges:
    let
      p = s.pos[e.a]
      q = s.pos[e.b]
      d = tangent(p, q)
      pz = dot(s.chart[2], p)
      dz = dot(s.chart[2], d)
      r = sqrt(pz*pz + dz*dz)
      u0 = arctan2(dz, pz)
    var cm = max(abs(pz), abs(dot(s.chart[2], q)))
    for k in -1..1:
      let u = u0 + PI*float(k)
      if u > 0.0 and u < e.len: cm = r
    result.link = min(result.link, arccos(min(1.0, cm)))
