import field

type
  RandomGeneratorBase* {.inheritable,pure.} = object
    count: int

func getCount*(x: var RandomGeneratorBase): int {.inline.} =
  result = x.count

func setCount*(x: var RandomGeneratorBase, c = 0) {.inline.} =
  x.count = c

func incCount*(x: var RandomGeneratorBase, c = 1) {.inline.} =
  x.count += c

type
  RandomGenerator* = concept var c
    c.type.high
    c.seed(1,1)
    c.next
    c.skip(1)

template high*[T:RandomGenerator](e: T): uint = (type T).high

func seed*[T:RandomGenerator](e: typedesc[T]; sed,index: auto): T =
  mixin seed
  result.seed(sed, index)

########################## MILC6 #############################
import milcrng

type
  Milc6Generator* = object of RandomGeneratorBase
    state*: RngMilc6

template high*(e: typedesc[Milc6Generator]): uint =
  mixin maxInt
  uint RngMilc6.maxInt

func seed*(e: var Milc6Generator; sed,index: auto) =
  e.setCount
  e.state.seedIndep(sed, index)

func milc6Generator*(sed,index: auto): Milc6Generator =
  result.seed(sed, index)

func next*(e: var Milc6Generator): uint =
  e.incCount
  result = uint e.state.integer

func skip*(e: var Milc6Generator, c = 1) =
  e.incCount(c)
  e.state.skip(c)

func seed*(x: Field[1, Milc6Generator]; sed: auto) =
  for e in x:
    x[e].seed(sed)

########################## MRG32k3a #############################
import mrg32k3a

type
  Mrg32k3aGenerator* = object of RandomGeneratorBase
    state*: MRG32k3a

template high*(e: typedesc[Mrg32k3aGenerator]): uint =
  uint MRG32k3a.maxInt

func seed*(e: var Mrg32k3aGenerator; sed,index: auto) =
  e.setCount
  e.state.seedIndep(sed, index)

func mrg32k3aGenerator*(sed,index: auto): Mrg32k3aGenerator =
  result.seed(sed, index)

func next*(e: var Mrg32k3aGenerator): uint =
  e.incCount
  result = uint e.state.integer

func skip*(e: var Mrg32k3aGenerator, c = 1) =
  e.incCount(c)
  e.state.skip(uint c)

########################## concat #############################
type
  ConcatGenerator*[T] = object
    state0*: T
template state*[T](e: ConcatGenerator[T]): T = e.state0
template state*[T](e: ConcatGenerator[ptr T]): T = e.state0[]
template count*[T](e: ConcatGenerator[T]): untyped = e.state.count
template `count=`*[T](e: var ConcatGenerator[T], c: auto) =
  e.state.count = c
template high*[T](e: typedesc[ConcatGenerator[T]]): uint = T.high*(T.high+2'u)
template high*[T](e: typedesc[ConcatGenerator[ptr T]]): uint = ConcatGenerator[T].high

func getCount*(x: var ConcatGenerator): int {.inline.} =
  result = x.count div 2

func setCount*(x: var ConcatGenerator, c = 0) {.inline.} =
  x.count = 2 * c

func incCount*(x: var ConcatGenerator, c = 1) {.inline.} =
  x.count += 2 * c

func seed*(e: var ConcatGenerator; sed,index: auto) =
  e.state.seed(sed, index)

func concatGenerator*[T:RandomGenerator](e: T): ConcatGenerator[T] =
  result.state0 = e

func concatGenerator*[T:RandomGenerator](e: ptr T): ConcatGenerator[ptr T] =
  result.state0 = e

func concatGenerator*[T](sed,index: auto): ConcatGenerator[T] =
  result.state.seed(sed, index)

func next*[T](e: var ConcatGenerator[T]): uint =
  const n = e.state.type.high + 1'u
  let v0 = e.state.next
  let v1 = e.state.next
  result = v0*n + v1

func skip*(e: var ConcatGenerator, c = 1) =
  e.state.skip(2*c)

################ tests #####################
when isMainModule:
  proc testE[T:RandomGenerator](e: var T) =
    echo e.getCount, " ", e.high
    echo e.next
    echo e.next
    echo e.getCount
    e.setcount
    echo e.getCount
    e.skip 5
    echo e.getCount
  proc test(T: typedesc) =
    var e = T.seed(987654321, 0)
    testE(e)
  proc testC(T: typedesc) =
    var e = T.seed(987654321, 0)
    testE(e)
    var c = concatGenerator(e)
    testE(c)
    testE(e)
    var v = concatGenerator(addr e)
    testE(v)
    testE(e)

  test(Milc6Generator)
  test(Mrg32k3aGenerator)
  test(ConcatGenerator[Milc6Generator])
  test(ConcatGenerator[Mrg32k3aGenerator])
  testC(Milc6Generator)
  testC(Mrg32k3aGenerator)
