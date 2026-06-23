suite "gauge coeffs":
  proc checkRectCoeff(coeff: Gactcoeff,
                      beta: float,
                      c1: float) =
    discard coeff.eval
    let c = coeff.cval
    check almostEqual(c.plaq, beta * (1.0 - 8.0 * c1))
    check almostEqual(c.rect, beta * c1)
    check almostEqual(c.pgm, 0.0)
    check almostEqual(c.adjplaq, 0.0)

  test "coefficient update refreshes cached coefficient expressions":
    let coeff = grt.toGvalue(GaugeActionCoeffs(plaq: 1.0, rect: 2.0))
    let coeffNorm = redot(coeff, coeff)

    coeffNorm :~ 5.0
    coeff.update GaugeActionCoeffs(plaq: 3.0, rect: 4.0)
    coeffNorm :~ 25.0

  test "direct coefficient field writes do not mark freshness":
    let coeff = grt.toGvalue(GaugeActionCoeffs(plaq: 1.0))
    let coeffNorm = redot(coeff, coeff)

    coeffNorm :~ 1.0
    coeff.cval = GaugeActionCoeffs(plaq: 3.0)
    check almostEqual(coeff.cval.plaq, 3.0)
    coeffNorm :~ 1.0

    coeff.update GaugeActionCoeffs(plaq: 3.0)
    coeffNorm :~ 9.0

  test "named action coefficient constructors preserve concrete type and values":
    let beta = grt.toGvalue(6.0)
    let wilson: Gactcoeff = actWilson(beta)
    let symanzik: Gactcoeff = actSymanzik(beta)
    let iwasaki: Gactcoeff = actIwasaki(beta)
    let dbw2: Gactcoeff = actDBW2(beta)
    let wilsonNorm: Gscalar = redot(wilson, wilson)

    discard wilson.eval
    check almostEqual(wilson.cval.plaq, 6.0)
    check almostEqual(wilson.cval.rect, 0.0)
    check almostEqual(wilson.cval.pgm, 0.0)
    check almostEqual(wilson.cval.adjplaq, 0.0)
    wilsonNorm :~ 36.0

    symanzik.checkRectCoeff(6.0, -1.0 / 12.0)
    iwasaki.checkRectCoeff(6.0, -0.331)
    dbw2.checkRectCoeff(6.0, -1.4088)

    beta.update 5.0
    symanzik.checkRectCoeff(5.0, -1.0 / 12.0)

  test "action coefficient erased copy compatibility stays direct":
    let source = grt.toGvalue(GaugeActionCoeffs(plaq: 2.0, rect: 3.0))
    let target = Gactcoeff(source.newOneOf)
    let scalarValue = grt.toGvalue(1.0)

    target.valCopy(source)
    check almostEqual(target.cval.plaq, source.cval.plaq)
    check almostEqual(target.cval.rect, source.cval.rect)
    check target.copyCompatible(source)
    check not target.copyCompatible(scalarValue)

  test "named action coefficient constant leaf starts fresh":
    let wilson: Gactcoeff = actWilson(scalar.toGvalue(grt, 6.0))
    let unitCoeff = wilson.inputs[1]

    check unitCoeff.epoch > 0

  test "actAdj keeps adjoint coefficient live":
    let beta = grt.toGvalue(6.0)
    let adjFac = grt.toGvalue(0.25)
    let coeff = actAdj(beta, adjFac)

    discard coeff.eval
    check almostEqual(coeff.cval.plaq, 6.0)
    check almostEqual(coeff.cval.rect, 0.0)
    check almostEqual(coeff.cval.pgm, 0.0)
    check almostEqual(coeff.cval.adjplaq, 1.5)

    adjFac.update 0.5
    discard coeff.eval
    check almostEqual(coeff.cval.plaq, 6.0)
    check almostEqual(coeff.cval.adjplaq, 3.0)

  test "actAdj scalarized gradients follow beta and adjoint factor":
    let beta = grt.toGvalue(6.0)
    let adjFac = grt.toGvalue(0.25)
    let seed = grt.toGvalue(GaugeActionCoeffs(plaq: 2.0, adjplaq: 3.0))
    let coeff = actAdj(beta, adjFac)
    let z = redot(coeff, seed)

    z :~ 16.5
    z.grad(beta) :~ 2.75
    z.grad(adjFac) :~ 18.0

    adjFac.update 0.5
    z :~ 21.0
    z.grad(beta) :~ 3.5

  test "coefficient backward hooks reject missing upstream":
    let beta = grt.toGvalue(6.0)
    let adjFac = grt.toGvalue(0.25)
    let coeff = actAdj(beta, adjFac)

    expect(GraphValueError):
      discard coeff.gfunc.backward(nil, coeff, 0, beta)
