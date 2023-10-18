import generator
import base/basicOps
import distributionUtils
export distributionUtils

########################## uniform int #####################################
type
  UniformInt* = object
    min: int
    max: int
func `$`*(d: UniformInt): string =
  result = $d.type & system.`$`(d)
func uniformInt*(min: int, max: int): UniformInt =
  result.min = min
  result.max = max
func generate*[T:RandomGenerator](d: UniformInt, g: var T): int =
  let n = uint(d.max - d.min + 1)
  let k = g.high div n
  let m = k * n
  var u = g.next
  while u >= m:
    u = g.next
  let t = u mod n
  result = d.min + int(t)

########################## uniform real #####################################
type
  UniformRealClosedClosed* = object
    min: float
    max: float
func `$`*(d: UniformRealClosedClosed): string =
  result = $d.type & system.`$`(d)
func uniformRealClosedClosed*(min: float, max: float): UniformRealClosedClosed =
  result.min = min
  result.max = max
func generate*[T:RandomGenerator](d: UniformRealClosedClosed, g: var T): float =
  const den = 1.0/float(g.high)
  let scale = (d.max-d.min)*den
  result = d.min + scale * float(g.next)

type
  UniformRealClosedOpen* = object
    min: float
    max: float
func `$`*(d: UniformRealClosedOpen): string =
  result = $d.type & system.`$`(d)
func uniformRealClosedOpen*(min: float, max: float): UniformRealClosedOpen =
  result.min = min
  result.max = max
func generate*[T:RandomGenerator](d: UniformRealClosedOpen, g: var T): float =
  const den = 1.0/(float(g.high)+1.0)
  let scale = (d.max-d.min)*den
  result = d.min + scale * g.next.float

type
  UniformRealOpenOpen* = object
    min: float
    max: float
func `$`*(d: UniformRealOpenOpen): string =
  result = $d.type & system.`$`(d)
func uniformRealOpenOpen*(min: float, max: float): UniformRealOpenOpen =
  result.min = min
  result.max = max
func generate*[T:RandomGenerator](d: UniformRealOpenOpen, g: var T): float =
  const den = 1.0/(float(g.high)+1.0)
  let scale = (d.max-d.min)*den
  result = d.min + scale * (float(g.next) + 0.5)

########################## normal #####################################

type
  NormalAlgorithm = enum
    naBoxMuller, naPolar, naCDF
  NormalDistribution* = object
    mean: float
    stddev: float
    alg: NormalAlgorithm
func `$`*(d: NormalDistribution): string =
  result = $d.type & system.`$`(d)
func normalDistribution*(mean = 0.0, stddev = 1.0, alg = naBoxMuller):
     NormalDistribution =
  result.mean = mean
  result.stddev = stddev
  result.alg = alg

proc gaussianBoxMuller*[T:RandomGenerator](g: var T): float =
  const h = float g.high
  const scale = -1.0/(h+1.0)
  const bias2 = 0.5*h
  const scale2 = PI/(h+1.0)
  let v = scale*(0.5 + float(g.next))
  let r = sqrt(-2.0 * ln1p(v))
  let p = scale2*(float(g.next) - bias2)
  result = r * sin(p)

proc erfcinv(y0: float): float =
  var y = y0
  if y0 > 1.0:
    y = 2.0 - y0
  var x = 0.0
  while true:
    let r = y - erfc(x)
    let d = 0.5*sqrt(PI)*exp(x*x)*r
    #echo "x: ", x, "  r: ", r, "  d: ", d
    x -= d
    if abs(x) > 26.0:
      break
    if abs(d) < 1e-14: break
  if y0 > 1.0:
    x = -x
  let r = erfc(x) - y0
  if abs(r) > 1e-15:
    echo "erfcinv: ", y0, "  ", x, "  ", r
  result = x

proc gaussianCDF*[T:RandomGenerator](g: var T): float =
  const h = float g.high
  const scale = 2.0/(h+1.0)
  let y = scale*(0.5 + float(g.next))
  result = sqrt(2.0)*erfcinv(y)

proc gaussianPolar*[T:RandomGenerator](g: var T): float =
  const h = float g.high
  const scale = 2.0/(h+1.0)
  const bias = 0.5*h
  while true:
    var x = scale*(g.next - bias)
    var y = scale*(g.next - bias)
    let rsq = x*x + y*y
    if rsq < 1.0 and rsq != 0.0:
      let fac = sqrt(-2.0 * ln(rsq) / rsq )
      result = x * fac
      break

func generate*[T:RandomGenerator](d: NormalDistribution, g: var T): float =
  var t = case d.alg
          of naBoxMuller: gaussianBoxMuller(g)
          of naPolar: gaussianPolar(g)
          of naCDF: gaussianBoxMuller(g)
  result = d.mean + d.stddev * t

######## z2, z4, ZN, U1, ... #####
  # Z_N  bias [0-1)   exp(i2pi*(k+bias)/N)    exp(i2pi*(k+bias-(h div 2))/N)

################ tests #####################
when isMainModule:
  proc test(d: auto, T: typedesc) =
    echo "****  ", d, "  ", T
    var e: T
    e.seed(987654321, 0)
    const n = 20
    var a: array[n,type(d.generate(e))]
    for i in 0..<n:
      a[i] = d.generate(e)
    echo a

  proc testD(d: auto) =
    test(d, Milc6Generator)
    test(d, Mrg32k3aGenerator)
    test(d, ConcatGenerator[Milc6Generator])
    test(d, ConcatGenerator[Mrg32k3aGenerator])

  proc testE(e: typedesc) =
    block:
      var d = uniformInt(-9, 9)
      test(d, e)
    block:
      var d = uniformRealClosedClosed(-1.0, 1.0)
      test(d, e)
    block:
      var d = uniformRealClosedOpen(-1.0, 1.0)
      test(d, e)
    block:
      var d = uniformRealOpenOpen(-1.0, 1.0)
      test(d, e)
    block:
      var d = normalDistribution(1.0, 2.0, naBoxMuller)
      test(d, e)
    block:
      var d = normalDistribution(1.0, 2.0, naCDF)
      test(d, e)
    block:
      var d = normalDistribution(1.0, 2.0, naPolar)
      test(d, e)

  testE(Milc6Generator)
  testE(ConcatGenerator[Milc6Generator])
  testE(Mrg32k3aGenerator)
  testE(ConcatGenerator[Mrg32k3aGenerator])
