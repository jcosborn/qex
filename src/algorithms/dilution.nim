import base, layout, field

type
  DilutionKind* = enum
    dkEvenOdd, dkCorners3D
  Dilution* = object
    case kind*: DilutionKind
    of dkEvenOdd: eo*: range[0..1]
    of dkCorners3D: c3d*: range[0..7]
template high(dk: DilutionKind): int =
  case dk
  of dkEvenOdd: 1
  of dkCorners3D: 7
template newDilution(dk: DilutionKind, i: int): Dilution =
  case dk
  of dkEvenOdd: Dilution(kind: dkEvenOdd, eo: i)
  of dkCorners3D: Dilution(kind: dkCorners3D, c3d: i)
proc `$`*(x: Dilution): string =
  case x.kind
  of dkEvenOdd: "EvenOdd " & $x.eo
  of dkCorners3D: "Corners3D " & $x.c3d

template sitesI(l: Layout, d: Dilution): auto =
  case d.kind
  of dkEvenOdd:
    # Assuming even-odd layout
    if d.eo == 0: itemsI(0, l.nEven)
    else: itemsI(l.nEven, l.nSites)
  of dkCorners3D:
    let
      n = l.nSites
      a = (threadNum*n) div numThreads
      b = (threadNum*n+n) div numThreads
      c = d.c3d
    var i = a
    while i < b:
      let
        x = l.coords[0][i].int and 1
        y = l.coords[1][i].int and 1
        z = l.coords[2][i].int and 1
      if (x + (y shl 1) + (z shl 2)) == c: yield i
      i.inc

iterator sites*(l: Layout, d: Dilution): int = l.sitesI d
iterator sites*(f: Field, d: Dilution): int = f.l.sitesI d

iterator dilution*(dl:DilutionKind): Dilution =
  #case dl
  #of dkEvenOdd:
  #  for i in 0..1:
  #    yield Dilution(kind:dkEvenOdd, eo:i)
  #of dkCorners3D:
  #  for i in 0..7:
  #    yield Dilution(kind:dkCorners3D, c3d:i)
  for i in 0..high(dl):
    yield newDilution(dl, i)

proc parseDilution*(dl:string): DilutionKind =
  case dl
  of "EO": return dkEvenOdd
  of "CORNER": return dkCorners3D
  else:
    echo "ERROR: unsupported dilution type: ",dl
    qexAbort()
