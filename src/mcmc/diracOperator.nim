import qex
import physics/[stagD, wilsonD]

import abstractFields
import typeUtilities

type
  DiracOperator*[S,T,U,V,W] = object
    case discretization: FieldType
      of StaggeredMatterField:
        stag*: Staggered[S,V]
        stagShifter*: seq[Shifter[T,V]]
        stagPsi*: T
      of WilsonMatterField:
        wils*: Wilson[S,W]
        wilsShifter*: seq[Shifter[U,W]]
        wilsPsi*: U
      else: discard

proc newDiracOperator[S,T,U,V,W](
    g: auto;
    discretization: FieldType;
    s: typedesc[S];
    t: typedesc[T];
    u: typedesc[U];
    v: typedesc[V];
    w: typedesc[W]
  ): DiracOperator[S,T,U,V,W] =
  let l = g[0].l
  result = DiracOperator[S,T,U,V,W](discretization: discretization)
  case result.discretization:
    of StaggeredMatterField:
      result.stagPsi = l.ColorVector()
      result.stag = newStag(g)
      result.stagShifter = newSeq[Shifter[T,V]](g[0].l.nDim)
      for mu in 0..<g.len: 
        result.stagShifter[mu] = newShifter(result.stagPsi, mu, 1)
    of WilsonMatterField: 
      result.wilsPsi = l.DiracFermion()
      result.wils = newWilson(g)
      #[
      result.wilsShifter = newSeq[Shifter[U,W]](g[0].l.nDim)
      for mu in 0..<g.len: 
        result.wilsShifter[mu] = newShifter(result.wilsPsi, mu, 1)
      ]#
    else: discard

proc newDiracOperator*(
    g: auto; 
    discretization: FieldType
  ): auto =
  
  let l = g[0].l

  result = g.newDiracOperator(
    discretization, 
    l.typeS, 
    l.typeT, 
    l.typeU, 
    l.typeV, 
    l.typeW
  )

if isMainModule:
  qexInit()
  var 
    lat = intSeqParam("lat", @[4, 4, 4, 4])
    lo = lat.newLayout(@[1, 1, 1, 1])
    u = lo.newGauge()
    sPhi = lo.ColorVector()
    wPhi = lo.DiracFermion()
  for discretization in @[StaggeredMatterField, WilsonMatterField]:
    var D = u.newDiracOperator(discretization)
    case discretization:
      of StaggeredMatterField: 
        for mu in 0..<u.len: discard D.stagShifter[mu] ^* sPhi
      of WilsonMatterField: discard 
      else: discard
  qexFinalize()