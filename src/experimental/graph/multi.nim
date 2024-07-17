import core, scalar

type
  Gint* {.final.} = ref object of Gvalue
    ival: int
  Gmulti* {.final.} = ref object of Gvalue
    mval: seq[Gvalue]

proc getint*(x: Gvalue): int = Gint(x).ival

proc `getint=`*(x: Gvalue, y: int) =
  let xs = Gint(x)
  xs.ival = y

proc update*(x: Gvalue, y: int) =
  x.getint = y
  x.updated

converter toGvalue*(x: int): Gvalue =
  result = Gint(ival: x)
  result.updated

method newOneOf*(x: Gint): Gvalue = Gint()
method valCopy*(z: Gint, x: Gint) = z.ival = x.ival

method `$`*(x: Gint): string = $x.ival

proc getmulti*(x: Gvalue): seq[Gvalue] = Gmulti(x).mval

proc `getmulti=`*(x: Gvalue, y: seq[Gvalue]) =
  let xs = Gmulti(x)
  xs.mval = y

proc update*(x: Gvalue, y: seq[Gvalue]) =
  x.getmulti = y
  x.updated

proc toGvalue*(x: seq[Gvalue]): Gmulti =
  # proc instead of converter to avoid converting seq
  result = Gmulti(mval: x)
  result.updated

method newOneOf*(x: Gmulti): Gvalue =
  let r = Gmulti(mval: newseq[Gvalue](x.mval.len))
  for i in 0..<x.mval.len:
    r.mval[i] = newOneOf x.mval[i]
  r

method valCopy*(z: Gmulti, x: Gmulti) =
  for i in 0..<z.mval.len:
    z.mval[i] = x.mval[i]

method `$`*(x: Gmulti): string = $x.mval

method `[]`*(x: Gvalue, i: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod($x & "[" & $i & "]")
method updateAt*(x: Gvalue, i: Gvalue, y: Gvalue): Gvalue {.base.} = raiseErrorBaseMethod("updateAt(" & $x & "," & $i & "," & $y & ")")

proc getAtmb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0:
    return z.inputs[0].updateAt(i, zb)
  else:
    raiseValueError("i must be 0, got: " & $i)

proc getAtmf(v: Gvalue) =
  let x = Gmulti(v.inputs[0])
  let i = v.inputs[1].getint
  v.valCopy x.mval[i]

let getAtm = newGfunc(forward = getAtmf, backward = getAtmb, name = "getAtm")

method `[]`*(x: Gmulti, i: Gint): Gvalue =
  result = newOneOf x.mval[i.ival]
  result.inputs = @[Gvalue(x), i]
  result.gfunc = getAtm

proc updateAtmb(zb: Gvalue, z: Gvalue, i: int, dep: Gvalue): Gvalue =
  case i
  of 0:
    return zb.updateAt(i, 0.0)  # TODO: need a zero method
  of 2:
    return zb[i]
  else:
    raiseValueError("i must be 0 or 2, got: " & $i)

proc updateAtmf(v: Gvalue) =
  let x = Gmulti(v.inputs[0])
  let i = v.inputs[1].getint
  let z = Gmulti(v)
  for k in 0..<z.mval.len:
    if k == i:
      z.mval[k] = v.inputs[2]
    else:
      z.mval[k] = x.mval[i]

let updateAtm = newGfunc(forward = updateAtmf, backward = updateAtmb, name = "updateAtm")

method updateAt*(x: Gmulti, i: Gint, y: Gvalue): Gvalue = Gmulti(mval: newseq[Gvalue](x.mval.len), inputs: @[Gvalue(x), i, y], gfunc: updateAtm)
