import qex
import layout

template typeS*(l: Layout): untyped =
  type(l.ColorMatrix())

template typeT*(l: Layout): untyped = 
  type(l.ColorVector())

template typeU*(l: Layout): untyped = 
  type(l.DiracFermion())

template typeV*(l: Layout): untyped =
  type(l.ColorVector()[0])

template typeW*(l: Layout): untyped =
  type(spproj1p(l.DiracFermion()[0]))