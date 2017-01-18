import qex
import qcdTypes
import rng
import math

qexInit()

var Lx = 64
var Lt = 64
var lat = [Lx,Lx,Lx,Lt]
var lo = newLayout(lat)

var vd = lo.ColorVector()
var vs = lo.ColorVectorS()
var vs2 = lo.ColorVectorS()

threads:
  for i in lo.sites:
    vd{i} := exp(-1e-2*i.float)
    vs{i} := exp(-1e-2*i.float)

  echo vd.norm2
  echo vs.norm2

  vs2 := vs + vd  # implicitly becomes vs2 := toSingle(toDouble(vs)+vd)
  echo vs2.norm2

  vs2 := toSingle(vs.toDouble + vd)
  echo vs2.norm2

  vs2 := vs + vd.toSingle  # explicitly single precision
  echo vs2.norm2

  vs2 := vs + 1.0  # implicitly becomes 1.0'f32
  echo vs2.norm2

qexFinalize()
