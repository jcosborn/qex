import base
import maths
import layout
import layout/[shifts]
import gauge/[gaugeAction, gaugeUtils]
import field/[fieldET]

proc redotP2(f1: SomeField; f2: SomeField2): auto =
  tic("redotP2[" & $f1.type & "," & $f2.type & "]")
  mixin redot, iredot, simdSum, items, toDouble, eval
  const nc = f1[0].nrows.float
  var d: evalType(toDouble(redot(f1[0],f2[0])))
  let t1 = f1
  let t2 = f2
  for x in items(t1): d += nc - redot(t1[x], t2[x])
  toc("local")
  result = simdSum(d)
  toc("simdSum")
  f1.l.threadRankSum(result)
  toc("threadRankSum")

template redot2*(f1: SomeAllField; f2: SomeAllField2): untyped =
  when declared(subsetObject): redotP2(f1[subsetObject], f2)
  elif declared(subsetString): redotP2(f1[subsetString], f2)
  else: redotP2(f1, f2)

proc pnorm2P(f: SomeField): auto =
  tic("pnorm2P[" & $type(f) & "]")
  mixin norm2, inorm2, simdSum, items, toDouble
  let offset = 2.0*f.l.nDim.float
  var n2: evalType(norm2(toDouble(f[0])))
  for x in items(f):
    inorm2(n2, toDouble(f[x]))
    n2 -= offset
  toc("local")
  result = simdSum(n2)
  toc("simdSum")
  f.l.threadRankSum(result)
  toc("threadRankSum")

template pnorm2*(f: SomeAllField): auto =
  when declared(subsetObject): pnorm2P(f[subsetObject])
  elif declared(subsetString): pnorm2P(f[subsetString])
  else: pnorm2P(f)

template pnorm2*(f: Subsetted): auto = pnorm2P(f)