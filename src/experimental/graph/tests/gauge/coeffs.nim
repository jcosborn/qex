suite "gauge coeffs":
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
