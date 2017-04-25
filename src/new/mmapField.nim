import fieldProxy
import ../stdUtils
import memfiles

type
  MmapFieldObj[T] = object
    n: int
    m: Memfile
    p: ptr cArray[T]
  MmapField[T] = FieldProxy[MmapFieldObj[T]]

proc newFieldImpl(x: MmapFieldObj, y: typedesc): MmapField[y] =
  discard
proc newMmapField[T](fn: string, n: int): MmapField[T] =
  result[].n = n
  let bytes = n*sizeof(T)
  result[].m = memfiles.open(fn, mode=fmReadWrite, newFileSize=bytes)
  result[].p = cast[ptr cArray[T]](result[].m.mem)
proc close[T](x: var MmapField[T]) =
  memfiles.close(x[].m)

template `[]`(x: MmapFieldObj, y: typed): untyped =
  x.p[y]
template `[]=`(x: MmapFieldObj, y: typed, z: untyped): untyped =
  x.p[y] = z
template `len`*(x: MmapFieldObj): untyped = x.n
#iterator indices(x: FieldObj): int =
#  countup(0,99)
#iterator indices(x: FieldObj): int =
#  countup(0,99)
template indices(x: MmapFieldObj): untyped = 0..<x.n
fieldScalarOverloads(SomeNumber)

when isMainModule:
  const n = 100
  var x = newMmapField[float]("x", n)
  var y = newMmapField[float]("y", n)
  var z = newMmapField[float]("z", n)
  template `:=`*(r: SomeNumber, x: int): untyped = r = (type(r))(x)
  template `+`*(r: SomeNumber, x: int): untyped = r + (type(r))(x)

  x := 1
  echo $x
  y := 2
  z := 3
  x := y + z + 10
  echo $x

  x.close
  y.close
  z.close
