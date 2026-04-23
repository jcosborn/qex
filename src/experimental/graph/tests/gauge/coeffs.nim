suite "gauge coeffs":
  proc checkRectCoeff(coeff: Gactcoeff,
                      beta: float,
                      c1: float) =
    discard coeff.eval
    let c = coeff.getactcoeff
    check almostEqual(c.plaq, beta * (1.0 - 8.0 * c1))
    check almostEqual(c.rect, beta * c1)
    check almostEqual(c.pgm, 0.0)
    check almostEqual(c.adjplaq, 0.0)

  test "named action coefficient constructors preserve concrete type and values":
    let beta = grt.toGvalue(6.0)
    let wilson: Gactcoeff = actWilson(beta)
    let symanzik: Gactcoeff = actSymanzik(beta)
    let iwasaki: Gactcoeff = actIwasaki(beta)
    let dbw2: Gactcoeff = actDBW2(beta)
    let direct: Gactcoeff = grt.actDBW2(6.0)
    let wilsonNorm: Gscalar = redot(wilson, wilson)

    discard wilson.eval
    check almostEqual(wilson.getactcoeff.plaq, 6.0)
    check almostEqual(wilson.getactcoeff.rect, 0.0)
    check almostEqual(wilson.getactcoeff.pgm, 0.0)
    check almostEqual(wilson.getactcoeff.adjplaq, 0.0)
    wilsonNorm :~ 36.0

    symanzik.checkRectCoeff(6.0, -1.0 / 12.0)
    iwasaki.checkRectCoeff(6.0, -0.331)
    dbw2.checkRectCoeff(6.0, -1.4088)
    direct.checkRectCoeff(6.0, -1.4088)

    beta.update 5.0
    symanzik.checkRectCoeff(5.0, -1.0 / 12.0)

  test "actAdj keeps adjoint coefficient live":
    let beta = grt.toGvalue(6.0)
    let adjFac = grt.toGvalue(0.25)
    let coeff = actAdj(beta, adjFac)

    discard coeff.eval
    check almostEqual(coeff.getactcoeff.plaq, 6.0)
    check almostEqual(coeff.getactcoeff.rect, 0.0)
    check almostEqual(coeff.getactcoeff.pgm, 0.0)
    check almostEqual(coeff.getactcoeff.adjplaq, 1.5)

    adjFac.update 0.5
    discard coeff.eval
    check almostEqual(coeff.getactcoeff.plaq, 6.0)
    check almostEqual(coeff.getactcoeff.adjplaq, 3.0)
