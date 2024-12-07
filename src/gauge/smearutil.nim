import base, layout

proc symStaple*(s: auto, alp: float, g1: auto, g2: auto,
                s1: auto, s2: auto, tm: auto, sm: auto) =
  # s = alp*[ shift from back mu(g1.adj * g2 * s1.field) + g1 * s2.field * s1.field.adj ]
  # g1: side link
  # g2: middle link
  # s1: g1 shifted from forward g2 direction
  # s2: g2 shifted from forward g1 direction
  # tm: temp gauge field
  # sm: shifter for backwards g1 direction
  tic()
  mixin adj
  tm := g1.adj * g2 * s1.field
  discard sm ^* tm
  s += alp * ( g1 * s2.field * s1.field.adj )
  s += alp * sm.field
  let nc = g1[0].nrows
  let siteFlops = float(nc*nc*((6*nc+2*(nc-1))*5+4*2))
  toc("symStaple", flops=siteFlops*g1.l.nSites)

proc symStapleDeriv*(f1, f2: auto;  # output
                     g1, g2: auto; s1, s2: auto;  # same as symStaple
                     c: auto, s: auto;  # chain and shift
                     tm1, tm2: auto;  # temporary fields
                     sm1, sm2: auto;  # shifts
                    ) =
  # f1: deriv wrt g1
  # f2: deriv wrt g2
  # s: c shifted from forward g1 direction
  # sm1: shifter from backwards g1 direction
  # sm2: shifter from backwards g2 direction
  tic()
  mixin adj
  # ∪ s.field.adj * g1.adj * g2 * s1.field
  # ∪† s.field * s1.field.adj * g2.adj * g1
  # ∩ c.adj * g1 * s2.field * s1.field.adj
  # ∩† c * s1.field * s2.field.adj * g1.adj
  tm1 := g1.adj * c * s1.field  # ∩†2  s2
  tm2 := g2.adj * g1 * s.field  # ∪†1  s1
  tm2 += c.adj * g1 * s2.field  # ∩3   s1
  discard sm1 ^* tm1
  discard sm2 ^* tm2
  f1 += g2 * s1.field * s.field.adj  # ∪1   g1
  f1 += c * s1.field * s2.field.adj  # ∩†3  g1
  f2 += g1 * s.field * s1.field.adj  # ∪†2  g2
  f2 += sm1.field
  f1 += sm2.field
  toc("symStapleDeriv")

