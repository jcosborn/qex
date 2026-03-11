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
  let siteFlops = float((4*(6*nc+2*(nc-1))+2*4)*nc*nc)
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
  threadBarrier()
  discard sm1 ^*! tm1
  discard sm2 ^*! tm2
  threadBarrier()
  f1 += g2 * s1.field * s.field.adj  # ∪1   g1
  f1 += c * s1.field * s2.field.adj  # ∩†3  g1
  f2 += g1 * s.field * s1.field.adj  # ∪†2  g2
  f2 += sm1.field
  f1 += sm2.field
  toc("symStapleDeriv")

template symStapleVJP*(f1, f2: auto;  # output
                       g1, g2: auto; s1, s2: auto;  # same as symStaple
                       c: auto, s: auto;  # chain and shift
                       tm1, tm2: auto;  # temporary fields
                       sm1, sm2: auto  # shifts
                      ) =
  ## Alias for symStapleDeriv using kernel naming.
  symStapleDeriv(f1, f2, g1, g2, s1, s2, c, s, tm1, tm2, sm1, sm2)

template shiftedAssign(sh, tmp, expr: untyped) =
  tmp := expr
  threadBarrier()
  discard sh ^*! tmp
  threadBarrier()

template shiftedAccumulate(dst, sh, tmp, expr: untyped) =
  shiftedAssign(sh, tmp, expr)
  dst += sh.field

proc symStapleHVP*(df1, df2, acc1, acc2: auto;  # output: directional second derivative wrt g1,g2,c; acc1,acc2 are pre-allocated temp fields
                      g1, g2: auto; s1, s2: auto;  # same as symStaple
                      c: auto, s: auto;            # chain and its shift
                      dg1, dg2: auto;              # input variation directions
                      dc: auto;                    # variation of the chain
                      ds1, ds2, dsc: auto;         # PRE-SHIFTED tangent shifters
                      tm1, tm2: auto;              # temporary fields
                      sm1, sm2: auto;              # shifts
                      eps: float = 0.0) =
  ## Analytical directional second derivative of the symmetric staple backprop.
  ##
  ## Given the first-derivative kernel
  ##   (f1,f2) = symStapleDeriv(g1,g2; c,s),
  ##
  ## this returns, in directions (dg1,dg2,dc),
  ##   (df1,df2) = d/dε [ symStapleDeriv(g1+ε dg1, g2+ε dg2; c+ε dc) ] |_{ε=0},
  ##
  ## including the dependence on the shifted fields s1 = shift(g1),
  ## s2 = shift(g2), and s = shift(c) using the same shifter geometry
  ## provided in the arguments.
  ##
  ## ds1, ds2, dsc are PRE-CREATED shifters for dg1, dg2, dc respectively.
  ## The caller must:
  ##   1. Create shifters OUTSIDE threads block with newShifter(field, dir, 1)
  ##   2. Apply shifts INSIDE threads block with ^*! before calling this function
  ##   3. Pass the shifter objects (not .field) to this function
  ##   4. Pass pre-allocated acc1, acc2 fields; no allocation inside threads.
  discard eps
  tic()
  mixin adj

  # Clear accumulators (pre-allocated by caller to avoid allocation inside threads)
  acc1 := 0
  acc2 := 0

  # Base shifted fields (inputs are expected to be pre-shifted).
  let s1b = s1.field
  let s2b = s2.field
  let sb  = s.field

  # Use passed shifters' .field for pre-shifted tangents.
  let ds1b = ds1.field
  let ds2b = ds2.field
  let dsb  = dsc.field

  # d(sm1) from dtm1 = dg1† c s1 + g1† c ds1 + g1† dc s1
  tm1 := dg1.adj * c
  tm1 := tm1 * s1b
  tm2 := g1.adj * c
  tm2 := tm2 * ds1b
  tm1 += tm2
  tm2 := g1.adj * dc
  tm2 := tm2 * s1b
  tm1 += tm2
  threadBarrier()
  discard sm1 ^*! tm1
  threadBarrier()

  # d(sm2) from dtm2 pieces:
  #   dg2† g1 s  +  g2† dg1 s  +  g2† g1 ds
  #   + dc† g1 s2 + c† dg1 s2 + c† g1 ds2
  tm1 := dg2.adj * g1
  tm1 := tm1 * sb
  tm2 := g2.adj * dg1
  tm2 := tm2 * sb
  tm1 += tm2
  tm2 := g2.adj * g1
  tm2 := tm2 * dsb
  tm1 += tm2
  tm2 := dc.adj * g1
  tm2 := tm2 * s2b
  tm1 += tm2
  tm2 := c.adj * dg1
  tm2 := tm2 * s2b
  tm1 += tm2
  tm2 := c.adj * g1
  tm2 := tm2 * ds2b
  tm1 += tm2
  threadBarrier()
  discard sm2 ^*! tm1
  threadBarrier()

  # acc1 contributions (was df1)
  tm1 := dg2 * s1b
  tm1 := tm1 * sb.adj
  acc1 += tm1

  tm1 := g2 * ds1b
  tm1 := tm1 * sb.adj
  acc1 += tm1

  tm1 := g2 * s1b
  tm1 := tm1 * dsb.adj
  acc1 += tm1

  tm1 := dc * s1b
  tm1 := tm1 * s2b.adj
  acc1 += tm1

  tm1 := c * ds1b
  tm1 := tm1 * s2b.adj
  acc1 += tm1

  tm1 := c * s1b
  tm1 := tm1 * ds2b.adj
  acc1 += tm1

  acc1 += sm2.field

  # acc2 contributions (was df2)
  tm1 := dg1 * sb
  tm1 := tm1 * s1b.adj
  acc2 += tm1

  tm1 := g1 * dsb
  tm1 := tm1 * s1b.adj
  acc2 += tm1

  tm1 := g1 * sb
  tm1 := tm1 * ds1b.adj
  acc2 += tm1

  acc2 += sm1.field

  # Add to outputs (no 0.5 factor - the derivative is exact)
  df1 += acc1
  df2 += acc2

  toc("symStapleHVP")

proc symStapleJVP*(ds: auto, alp: float,
                       g1, g2: auto, s1, s2: auto,
                       dg1, dg2: auto, ds1, ds2: auto,
                       tm: auto, sm: auto) =
  ## Forward tangent of symStaple.
  ##
  ## Computes: ds += alp * d/dε[ symStaple(g1+ε dg1, g2+ε dg2, ...) ]|_{ε=0}
  ##
  ## Recall symStaple computes:
  ##   s += alp * [ shift_back(g1† g2 s1) + g1 s2 s1† ]
  ##
  ## Tangent:
  ##   ds += alp * [ shift_back(dg1† g2 s1 + g1† dg2 s1 + g1† g2 ds1)
  ##               + dg1 s2 s1† + g1 ds2 s1† + g1 s2 ds1† ]
  tic()
  mixin adj
  # Term 1: shift_back(dg1† g2 s1 + g1† dg2 s1 + g1† g2 ds1)
  tm := dg1.adj * g2 * s1.field
  tm += g1.adj * dg2 * s1.field
  tm += g1.adj * g2 * ds1.field
  discard sm ^* tm
  ds += alp * sm.field

  # Term 2: dg1 s2 s1† + g1 ds2 s1† + g1 s2 ds1†
  ds += alp * (dg1 * s2.field * s1.field.adj)
  ds += alp * (g1 * ds2.field * s1.field.adj)
  ds += alp * (g1 * s2.field * ds1.field.adj)
  toc("symStapleJVP")

proc symStapleJVPVJP*(dg1bar, dg2bar: auto, alp: float,
                          g1, g2: auto, s1, s2: auto,
                          dsbar: auto, shiftedDsbar: auto,
                          tm1, tm2: auto,
                          sm1, sm2: auto) =
  ## Adjoint of symStapleJVP w.r.t. (dg1, dg2).
  ##
  ## Given dsbar (the adjoint seed), accumulates into (dg1bar, dg2bar) such that:
  ##   <dsbar, symStapleJVP(dg1, dg2, ds1, ds2, ...)> = <dg1bar, dg1> + <dg2bar, dg2>
  ##
  ## This is the COMPLETE adjoint including terms from shifted tangents ds1, ds2.
  ##
  ## Parameters:
  ##   dg1bar, dg2bar: output (accumulated)
  ##   g1, g2: gauge links (same as symStapleJVP)
  ##   s1, s2: PRE-SHIFTED gauge link shifters, must be freshly applied before calling:
  ##           discard s1 ^*! g1; discard s2 ^*! g2; threadBarrier()
  ##   dsbar: adjoint seed
  ##   shiftedDsbar: scratch forward shifter in the g1 direction, reused to shift dsbar
  ##   tm1, tm2: temporary fields
  ##   sm1: backward shifter in g1 direction (same as symStapleJVP's sm)
  ##   sm2: backward shifter in g2 direction
  ##
  ## Recall symStapleJVP computes:
  ##   ds += alp * [ shift_back_ν(dg1† g2 s1 + g1† dg2 s1 + g1† g2 ds1)
  ##               + dg1 s2 s1† + g1 ds2 s1† + g1 s2 ds1† ]
  ##
  ## where ν is the g1 direction, μ is the g2 direction,
  ## s1 = g1 shifted forward in μ, s2 = g2 shifted forward in ν,
  ## ds1 = dg1 shifted forward in μ, ds2 = dg2 shifted forward in ν.
  ##
  ## The adjoint has the SAME STRUCTURE as symStapleDeriv but with different field roles:
  ## - In symStapleDeriv: chain c is the linear input
  ## - Here: (dg1, dg2) are the linear inputs, and we treat dsbar as the chain
  ##
  ## The complete adjoint is:
  ##   dg1bar += alp * [ g2 s1 (shift_fwd_ν dsbar)† + dsbar s1 s2† ]
  ##           + alp * shift_back_μ[ g2 g1 (shift_fwd_ν dsbar) + dsbar† g1 s2 ]
  ##   dg2bar += alp * g1 (shift_fwd_ν dsbar) s1†
  ##           + alp * shift_back_ν[ g1 dsbar s1† ]
  tic()
  mixin adj

  shiftedAssign(shiftedDsbar, tm1, dsbar)

  # === Terms from dg1 directly (A, B) ===
  # Term A: ds[x] += alp * dg1[x+ν]† g2[x+ν] s1[x+ν]  (from shift_back in ν)
  # Shift dsbar forward in ν once, then apply the local adjoint
  # g2 * s1 * shift_fwd_ν(dsbar)†.
  tm2 := g2 * s1.field * shiftedDsbar.field.adj
  dg1bar += alp * tm2

  # Term B: ds[x] += alp * dg1[x] s2_orig[x] s1[x]†
  # Adjoint: dg1bar[x] += alp * dsbar[x] s1.field[x] s2_orig.field[x]†
  #        = alp * dsbar[x] s1.field[x] s2.field[x]†
  tm2 := dsbar * s1.field * s2.field.adj
  dg1bar += alp * tm2

  # === Terms from dg2 directly (C) ===
  # Term C: ds[x] += alp * g1[x+ν]† dg2[x+ν] s1[x+ν]  (from shift_back in ν)
  # Adjoint: dg2bar[y] += alp * g1[y] dsbar[y-ν] s1.field[y]†
  #        = alp * g1[y] shiftedDsbar.field[y] s1.field[y]†
  tm2 := g1 * shiftedDsbar.field * s1.field.adj
  dg2bar += alp * tm2

  # === Terms from ds1 = shifted dg1 (D, F) ===
  # These two terms contribute backward shifts in μ.
  #
  # Term D: shift_back_μ(g2† * g1 * shift_fwd_ν(dsbar))
  # Then shift backward in μ using sm2.
  shiftedAssign(sm2, tm2, g2.adj * g1 * shiftedDsbar.field)  # sm2.field[z] = tm2[z+μ]
  dg1bar += alp * sm2.field

  # Term F: shift_back_μ(dsbar† * g1 * s2)
  # Then shift backward in μ.
  shiftedAssign(sm2, tm2, dsbar.adj * g1 * s2.field)
  dg1bar += alp * sm2.field

  # === Term from ds2 = shifted dg2 (E) ===
  # Term E: shift_back_ν(g1† * dsbar * s1)
  # Then shift backward in ν using sm1.
  shiftedAssign(sm1, tm2, g1.adj * dsbar * s1.field)
  dg2bar += alp * sm1.field

  toc("symStapleJVPVJP")

proc symStapleHVP*(df1, df2, acc1, acc2: auto;  # output + pre-allocated accumulators
                      g1, g2: auto; s1, s2: auto;  # same as symStaple
                      c: auto, s: auto;            # chain and shift (held fixed)
                      dg1, dg2: auto;              # input variation directions
                      ds1, ds2, dsc: auto;         # PRE-SHIFTED tangent shifters
                      tm1, tm2: auto;              # temporary fields
                      sm1, sm2: auto;              # shifts
                      eps: float = 0.0) =
  ## Compatibility wrapper when no upstream chain variation dc is provided.
  ## Creates a zero dc field internally.
  var dc {.noInit.}: type(c)
  dc := 0
  symStapleHVP(df1, df2, acc1, acc2, g1, g2, s1, s2, c, s, dg1, dg2, dc, ds1, ds2, dsc, tm1, tm2, sm1, sm2, eps)

proc symStapleVJPChain*(cbar: auto;  # output: adjoint wrt chain (accumulated)
                               g1, g2: auto; s1, s2: auto;  # geometry (same as symStapleDeriv)
                               f1bar, f2bar: auto;  # adjoint seeds
                               shiftedF2bar, shiftedF1bar: auto;  # scratch forward shifters in g1/g2 directions
                               tm1, tm2: auto;  # temporaries
                               sm1: auto  # backward shifter in g1 direction
                              ) =
  ## Adjoint of symStapleDeriv with respect to the chain field c.
  ##
  ## Given (f1bar, f2bar), computes cbar such that:
  ##   ⟨(f1bar,f2bar), symStapleDeriv(g1,g2,c)⟩ = ⟨cbar, c⟩
  ##
  ## s1: forward shifter in g2=μ direction (initially holds g1 shifted)
  ## s2: forward shifter in g1=ν direction (initially holds g2 shifted)
  ## shiftedF2bar shifts `f2bar` in the g1 direction, shiftedF1bar shifts `f1bar` in the g2 direction.
  tic()
  mixin adj

  # s1.field[x] = g1[x-μ] where μ = g2 direction (forward shift gets value from behind)
  # s2.field[x] = g2[x-ν] where ν = g1 direction

  # === Term 1 (line 47): f1[x] += c[x] s1[x] s2[x]† ===
  # cbar[x] += f1bar[x] s2[x] s1[x]†
  tm1 := f1bar * s2.field * s1.field.adj
  cbar += tm1

  # === Term 2 (line 46): f1[x] += g2[x] s1g[x] c[x-ν]† ===
  # For f = M c†, adjoint is cbar = fbar† M
  # cbar[z] = f1bar[z+ν]† g2[z+ν] g1[z+ν-μ]
  # Compute expr[w] = f1bar[w]† g2[w] s1.field[w] at site w
  # Then cbar[z] = expr[z+ν] → BACKWARD shift in ν (sm1)
  shiftedAccumulate(cbar, sm1, tm1, f1bar.adj * g2 * s1.field)

  # === Term 3 (line 48): f2[x] += g1[x] c[x-ν] s1g[x]† ===
  # Adjoint: cbar[z] = g1[z+ν]† f2bar[z+ν] g1[z+ν-μ]
  # Compute expr[w] = g1[w]† f2bar[w] s1.field[w] at site w
  # Then cbar[z] = expr[z+ν] → BACKWARD shift in ν (sm1)
  shiftedAccumulate(cbar, sm1, tm1, g1.adj * f2bar * s1.field)

  # === Term 6 (line 49 from 39): f2[x] += g1[x+ν]† c[x+ν] g1[x+ν-μ] ===
  # Adjoint: cbar[z] = g1[z] f2bar[z-ν] g1[z-μ]†
  shiftedAssign(shiftedF2bar, tm2, f2bar)
  tm1 := g1 * shiftedF2bar.field * s1.field.adj
  cbar += tm1

  # === Term 4 (line 50A from 40): f1[x] += g2[x+μ]† g1[x+μ] c[x+μ-ν] ===
  # Adjoint: cbar[z] = g1[z+ν]† g2[z+ν] f1bar[z-μ+ν]
  shiftedAssign(shiftedF1bar, tm2, f1bar)
  shiftedAccumulate(cbar, sm1, tm1, g1.adj * g2 * shiftedF1bar.field)

  # === Term 5 (line 50B from 41): f1[x] += c[x+μ]† g1[x+μ] g2[x+μ-ν] ===
  # For f = c† M, adjoint is cbar = M fbar†
  # cbar[z] = g1[z] g2[z-ν] f1bar[z-μ]†
  tm1 := g1 * s2.field * shiftedF1bar.field.adj
  cbar += tm1

  toc("symStapleVJPChain")

proc symStapleHVPVJP*(dg1bar, dg2bar, dcbar: auto;  # output (accumulated)
                          g1, g2: auto; s1, s2: auto;   # gauge links (fixed)
                          c: auto, s: auto;              # chain (fixed)
                          f1bar, f2bar: auto;            # adjoint seeds
                          shiftedF2bar, shiftedF1bar: auto;  # scratch forward shifters in g1/g2 directions
                          tm1, tm2: auto;                # temps
                          sm1, sm2: auto) =              # shifters
  ## Adjoint of symStapleHVP w.r.t. (dg1, dg2, dc).
  ##
  ## Given (f1bar, f2bar), accumulates into (dg1bar, dg2bar, dcbar) such that:
  ##   <(f1bar, f2bar), symStapleHVP(dg1, dg2, dc)> = <dg1bar, dg1> + <dg2bar, dg2> + <dcbar, dc>
  ##
  ## s1, s2, s: PRE-SHIFTED shifters (same as symStapleHVP)
  ##   s1.field = g1 shifted, s2.field = g2 shifted, s.field = c shifted
  ## sm1: backward shifter in g1 direction
  ## sm2: backward shifter in g2 direction
  ##
  ## Recall symStapleHVP computes (from lines 98-175 in the forward pass):
  ##   Part 1: sm1 terms from dg1†*c*s1 + g1†*c*ds1 + g1†*dc*s1
  ##   Part 2: sm2 terms from dg2†*g1*s + g2†*dg1*s + g2†*g1*ds + dc†*g1*s2 + c†*dg1*s2 + c†*g1*ds2
  ##   Part 3: acc1 terms (local)
  ##   Part 4: acc2 terms (local)
  ##
  ## The adjoint reverses all these operations.
  tic()
  mixin adj

  # Read pre-shifted fields (these are valid at entry)
  # s1.field = g1 shifted in μ (g2 direction)
  # s2.field = g2 shifted in ν (g1 direction)
  # s.field = c shifted in ν (g1 direction)

  # === Adjoint of acc1 local terms (lines 136-158) ===

  # Term 1 (line 136-138): acc1 += dg2 * s1b * sb†
  # Linear in dg2: adjoint is dg2bar += f1bar * sb * s1b†
  tm1 := f1bar * s.field * s1.field.adj
  dg2bar += tm1

  # Term 2 (line 140-142): acc1 += g2 * ds1b * sb†
  # Linear in ds1b = shifted dg1: adjoint is dg1bar += shift_back_μ(g2† * f1bar * sb)
  tm1 := g2.adj * f1bar * s.field
  threadBarrier()
  discard sm2 ^*! tm1  # sm2 shifts back in g2 (μ) direction
  threadBarrier()
  dg1bar += sm2.field

  # Term 3 (line 144-146): acc1 += g2 * s1b * dsb†
  # Linear in dsb† where dsb = shifted dc
  # For f = M * A†: <fbar, M A†> = <fbar† M, A>
  # M = g2 * s1b, A = dc shifted in ν
  # Adjoint: dcbar[z] = (fbar† M)[z+ν] = shift_back_ν(f1bar† * g2 * s1.field)
  tm1 := f1bar.adj * g2 * s1.field
  threadBarrier()
  discard sm1 ^*! tm1
  threadBarrier()
  dcbar += sm1.field

  # Term 4 (line 148-150): acc1 += dc * s1b * s2b†
  # Linear in dc: adjoint is dcbar += f1bar * s2b * s1b†
  tm1 := f1bar * s2.field * s1.field.adj
  dcbar += tm1

  # Term 5 (line 152-154): acc1 += c * ds1b * s2b†
  # Linear in ds1b = shifted dg1: adjoint is dg1bar += shift_back_μ(c† * f1bar * s2b)
  tm1 := c.adj * f1bar * s2.field
  threadBarrier()
  discard sm2 ^*! tm1
  threadBarrier()
  dg1bar += sm2.field

  # Term 6 (line 156-158): acc1 += c * s1b * ds2b†
  # Linear in ds2b† where ds2b = shifted dg2
  # For f = M * A†: <fbar, M A†> = <fbar† M, A>
  # M = c * s1b, A = dg2 shifted in ν
  # Adjoint: dg2bar[z] = (fbar† M)[z+ν] = shift_back_ν(f1bar† * c * s1.field)
  tm1 := f1bar.adj * c * s1.field
  threadBarrier()
  discard sm1 ^*! tm1
  threadBarrier()
  dg2bar += sm1.field

  # === Adjoint of acc2 local terms (lines 163-173) ===

  # Term 8 (line 163-165): acc2 += dg1 * sb * s1b†
  # Linear in dg1: adjoint is dg1bar += f2bar * s1b * sb†
  tm1 := f2bar * s1.field * s.field.adj
  dg1bar += tm1

  # Term 9 (line 167-169): acc2 += g1 * dsb * s1b†
  # Linear in dsb = shifted dc: adjoint is dcbar += shift_back_ν(g1† * f2bar * s1b)
  tm1 := g1.adj * f2bar * s1.field
  threadBarrier()
  discard sm1 ^*! tm1
  threadBarrier()
  dcbar += sm1.field

  # Term 10 (line 171-173): acc2 += g1 * sb * ds1b†
  # Linear in ds1b† where ds1b = shifted dg1
  # For f = M * A†: <fbar, M A†> = <fbar† M, A>
  # M = g1 * sb, A = dg1 shifted in μ
  # Adjoint: dg1bar[z] = (fbar† M)[z+μ] = shift_back_μ(f2bar† * g1 * s.field)
  tm1 := f2bar.adj * g1 * s.field
  threadBarrier()
  discard sm2 ^*! tm1
  threadBarrier()
  dg1bar += sm2.field

  # === Adjoint of sm1 shifted term (lines 98-109, contributes to acc2 via line 175) ===
  # Forward: tm1 = dg1† * c * s1b + g1† * c * ds1b + g1† * dc * s1b
  # sm1.field = shift_back_ν(tm1), then acc2 += sm1.field
  #
  shiftedAssign(shiftedF2bar, tm2, f2bar)

  # For term dg1† * c * s1b shifted by sm1:
  # <f2bar, shift_back_ν(dg1† * c * s1b)> = <shift_fwd_ν(f2bar), dg1† * c * s1b>
  # For A†: <sf, A† M> = <M sf†, A>
  # where M = c * s1b, sf = shift_fwd_ν(f2bar) = shiftedF2bar.field, A = dg1
  # dg1bar += M * sf† = c * s1.field * shiftedF2bar.field†
  tm1 := c * s1.field * shiftedF2bar.field.adj
  dg1bar += tm1

  # For term g1† * c * ds1b in sm1:
  # <shift_fwd_ν(f2bar), g1† * c * ds1b> where ds1b = dg1 shifted in μ
  # = <c† g1 sf, ds1b> where sf = shiftedF2bar.field
  # = <shift_back_μ(c† g1 sf), dg1>
  tm1 := c.adj * g1 * shiftedF2bar.field
  threadBarrier()
  discard sm2 ^*! tm1
  threadBarrier()
  dg1bar += sm2.field

  # For term g1† * dc * s1b in sm1:
  # <shift_fwd_ν(f2bar), g1† * dc * s1b> where sf = shiftedF2bar.field
  # Linear in dc (middle position): L(dc) = g1† * dc * s1b
  # L*(Y) = (g1†)† * Y * (s1b)† = g1 * Y * s1b† = g1 * sf * s1.field†
  # dcbar += g1 * shiftedF2bar.field * s1.field†
  tm1 := g1 * shiftedF2bar.field * s1.field.adj
  dcbar += tm1

  # === Adjoint of sm2 shifted term (lines 111-133, contributes to acc1 via line 160) ===
  # Forward: tm1 = dg2† * g1 * sb + g2† * dg1 * sb + g2† * g1 * dsb
  #              + dc† * g1 * s2b + c† * dg1 * s2b + c† * g1 * ds2b
  # sm2.field = shift_back_μ(tm1), then acc1 += sm2.field
  #
  shiftedAssign(shiftedF1bar, tm2, f1bar)

  # s.field still has c shifted (sb)
  # s2.field still holds the shifted g2 links from the caller.

  # For dg2† * g1 * sb:
  # <f1bar, shift_back_μ(dg2† * g1 * sb)> = <shift_fwd_μ(f1bar), dg2† * g1 * sb>
  # For A†: <sf, A† M> = <M sf†, A>
  # where M = g1 * sb, sf = shiftedF1bar.field, A = dg2
  # dg2bar += M * sf† = g1 * s.field * shiftedF1bar.field†
  tm1 := g1 * s.field * shiftedF1bar.field.adj
  dg2bar += tm1

  # For g2† * dg1 * sb:
  # <sf, g2† * dg1 * sb> where sf = shiftedF1bar.field (f1bar shifted forward in μ)
  # Linear in dg1 (middle position): L(dg1) = g2† * dg1 * sb
  # L*(Y) = (g2†)† * Y * (sb)† = g2 * Y * sb† = g2 * sf * s.field†
  # dg1bar += g2 * shiftedF1bar.field * s.field†
  tm1 := g2 * shiftedF1bar.field * s.field.adj
  dg1bar += tm1

  # For g2† * g1 * dsb where dsb = dc shifted in ν:
  # <sf, g2† * g1 * dsb> = <g1† g2 sf, dsb> = <shift_back_ν(g1† g2 sf), dc>
  tm1 := g1.adj * g2 * shiftedF1bar.field
  threadBarrier()
  discard sm1 ^*! tm1
  threadBarrier()
  dcbar += sm1.field

  # For dc† * g1 * s2b:
  # <sf, dc† * g1 * s2b> where sf = shiftedF1bar.field, s2b = s2.field = g2 shifted
  # For A†: <sf, A† M> = <M sf†, A>
  # where M = g1 * s2b, sf = shiftedF1bar.field, A = dc
  # dcbar += M * sf† = g1 * s2.field * shiftedF1bar.field†
  tm1 := g1 * s2.field * shiftedF1bar.field.adj
  dcbar += tm1

  # For c† * dg1 * s2b:
  # <sf, c† * dg1 * s2b> where sf = shiftedF1bar.field
  # Linear in dg1 (middle position): L(dg1) = c† * dg1 * s2b
  # L*(Y) = (c†)† * Y * (s2b)† = c * Y * s2b† = c * sf * s2.field†
  # dg1bar += c * shiftedF1bar.field * s2.field†
  tm1 := c * shiftedF1bar.field * s2.field.adj
  dg1bar += tm1

  # For c† * g1 * ds2b where ds2b = dg2 shifted in ν:
  # <sf, c† * g1 * ds2b> = <g1† c sf, ds2b> = <shift_back_ν(g1† c sf), dg2>
  tm1 := g1.adj * c * shiftedF1bar.field
  threadBarrier()
  discard sm1 ^*! tm1
  threadBarrier()
  dg2bar += sm1.field

  toc("symStapleHVPVJP")
