import base
import matrixConcept
import types
import matinv
import macros, algorithm

## polynomial power series exponential
## optimized for minimal number of matrix multiplies

template C2 = 0.5
template C3 = 1.0/6.0
template C4 = 1.0/24.0
template C5 = 1.0/120.0
template C6 = 1.0/720.0
template C7 = 1.0/5040.0
template C8 = 1.0/40320.0
template C9 = 1.0/362880.0
template C10 = 1.0/3628800.0
template C11 = 1.0/39916800.0
template C12 = 1.0/479001600.0

macro splitVars(sl: untyped): untyped =
  #echo "====="
  #echo sl.repr
  #echo "-----"
  result = newStmtList()
  sl.expectKind(nnkStmtList)
  for vs in sl:
    vs.expectKind(nnkVarSection)
    for id in vs:
      #echo id.repr
      id.expectKind(nnkIdentDefs)
      let x = id[0]
      id[1].expectKind(nnkEmpty)
      var e = id[2]
      var s = newSeq[NimNode]()
      while e.kind == nnkInfix and $e[0]=="+":
        let i = infix(x, "+=", e[2])
        #echo i.repr
        s.add i
        e = e[1]
      result.add newVarStmt(x, e)
      s.reverse
      result.add s
      #echo id.treerepr
    #var t = newVarStmt()
  #echo result.repr
  #return sl

# deriv convention: (d/dm') Tr[ w' f(m) + f(m)' w ]/2

# 2 MM
proc expm1Poly3*(m: Mat1): auto {.noInit.} =
  splitVars:
    var b = C3*m + C2
    var a =  b*m + 1
    var r =  a*m
  r
proc expPoly3*(m: Mat1): auto {.noInit.} =
  var r = expm1Poly3(m)
  r += 1
  r

# (d/dm') Tr[ w'(1+m+C2m^2+C3m^3) + (1+m+C2m^2+C3m^3)'w ] /2
# = w + C2(m'w+wm') + C3(m'm'w+m'wm'+wm'm')
# = w(1+m'(C2+C3m')) + m'(w(C2+C3m') + m'(C3w))
proc expm1Poly3Deriv*(m: Mat1, w: Mat2): auto {.noInit.} =
  splitVars:
    var a = C3*m.adj
    var b = a + C2        # b = C3m' + C2
    var d = w*b           # d = C3wm' + C2w
    var c = d*m.adj + w   # c = C3wm'm'+C2wm'+w
    var h = a*w + d       # h = C3(m'w+wm')+C2w
    var r = m.adj*h + c
  r
proc expPoly3Deriv*(m: Mat1, w: Mat2): auto {.noInit.} =
  expm1Poly3Deriv(m, w)

# 2 MM
proc expm1Poly4*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  splitVars:
    var a = C4*m2 + C3*m + C2
    var r =  a*m2 + m
  r
proc expPoly4*(m: Mat1): auto {.noInit.} =
  var r = expm1Poly4(m)
  r += 1
  r

# m2 = m' m'
# a = w m'
# b = m' w
# c = m' w m'
# d = wm'+m'w
# e = C4 d + C3 w
# f = m2 e + e m2
# r = C2 d + C3 c + f
# w + C2(m'w+wm') + C3(m'm'w+m'wm'+wm'm') + C4(m'm'm'w+m'm'wm'+m'wm'm'+wm'm'm')
proc expm1Poly4Deriv*(m: Mat1, w: Mat2): auto {.noInit.} =
  splitVars:
    #var a = C3*m.adj
    #var b = a + C2        # b = C3m' + C2
    var g = C4*m.adj
    var f = g + C3
    var e = w*f
    var a = g*w + e
    var d = e*m.adj + C2*w
    var c = d*m.adj + w
    var h = m.adj*a + d
    var r = m.adj*h + c
    # h = C4(m'm'w+m'wm'+wm'm') + C3(m'w+wm') + C2w
    # c = C4wm'm'm' + C3wm'm' + C2wm' + w
    # d = C4wm'm' + C3wm' + C2w
    # e = C4wm' + C3w
    # a = C4(m'w+wm') + C3w
  r
proc expPoly4Deriv*(m: Mat1, w: Mat2): auto {.noInit.} =
  expm1Poly4Deriv(m, w)

# 3 MM
proc expm1Poly5*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  splitVars:
    var b =        C5*m + C4
    var a = b*m2 + C3*m + C2
    var r = a*m2 + m
  r
# 3 MM
proc expPoly5*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  splitVars:
    var b =        C5*m + C4
    var a = b*m2 + C3*m + C2
    var r = a*m2 + m + 1
  r

# 3 MM
proc expm1Poly6*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  splitVars:
    var b = C6*m2 + C5*m + C4
    var a =  b*m2 + C3*m + C2
    var r =  a*m2 + m
  r
# 3 MM
proc expPoly6*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  splitVars:
    var b = C6*m2 + C5*m + C4
    var a =  b*m2 + C3*m + C2
    var r =  a*m2 + m + 1
  r

# 4 MM
proc expm1Poly7*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  splitVars:
    var c =        C7*m + C6
    var b = c*m2 + C5*m + C4
    var a = b*m2 + C3*m + C2
    var r = a*m2 + m
  r
# 4 MM
proc expPoly7*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  splitVars:
    var c =        C7*m + C6
    var b = c*m2 + C5*m + C4
    var a = b*m2 + C3*m + C2
    var r = a*m2 + m + 1
  r

# 4 MM
proc expm1Poly8*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  splitVars:
    var c = C8*m2 + C7*m + C6
    var b =  c*m2 + C5*m + C4
    var a =  b*m2 + C3*m + C2
    var r =  a*m2 + m
  r
# 4 MM
proc expPoly8*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  splitVars:
    var c = C8*m2 + C7*m + C6
    var b =  c*m2 + C5*m + C4
    var a =  b*m2 + C3*m + C2
    var r =  a*m2 + m + 1
  r

# 4 MM
proc expm1Poly9*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  let m3 = m2*m
  splitVars:
    var b = C9*m3 + C8*m2 + C7*m + C6
    var a =  b*m3 + C5*m2 + C4*m + C3
    var r =  a*m3 + C2*m2 + m
  r
# 4 MM
proc expPoly9*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  let m3 = m2*m
  splitVars:
    var b = C9*m3 + C8*m2 + C7*m + C6
    var a =  b*m3 + C5*m2 + C4*m + C3
    var r =  a*m3 + C2*m2 + m + 1
  r

# 5 MM
proc expm1Poly10*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  let m3 = m2*m
  splitVars:
    var c =               C10*m + C9
    var b = c*m3 + C8*m2 + C7*m + C6
    var a = b*m3 + C5*m2 + C4*m + C3
    var r = a*m3 + C2*m2 + m
  r
# 5 MM
proc expPoly10*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  let m3 = m2*m
  splitVars:
    var c =               C10*m + C9
    var b = c*m3 + C8*m2 + C7*m + C6
    var a = b*m3 + C5*m2 + C4*m + C3
    var r = a*m3 + C2*m2 + m + 1
  r

# 5 MM
proc expm1Poly11*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  let m3 = m2*m
  splitVars:
    var c =       C11*m2 + C10*m + C9
    var b = c*m3 + C8*m2 + C7*m + C6
    var a = b*m3 + C5*m2 + C4*m + C3
    var r = a*m3 + C2*m2 + m
  r
# 5 MM
proc expPoly11*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  let m3 = m2*m
  splitVars:
    var c =       C11*m2 + C10*m + C9
    var b = c*m3 + C8*m2 + C7*m + C6
    var a = b*m3 + C5*m2 + C4*m + C3
    var r = a*m3 + C2*m2 + m + 1
  r

# 5 MM
proc expm1Poly12*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  let m3 = m2*m
  splitVars:
    var c = C12*m3 + C11*m2 + C10*m + C9
    var b =   c*m3 + C8*m2 + C7*m + C6
    var a =   b*m3 + C5*m2 + C4*m + C3
    var r =   a*m3 + C2*m2 + m
  r
# 5 MM
proc expPoly12*(m: Mat1): auto {.noInit.} =
  let m2 = m*m
  let m3 = m2*m
  splitVars:
    var c = C12*m3 + C11*m2 + C10*m + C9
    var b =   c*m3 + C8*m2 + C7*m + C6
    var a =   b*m3 + C5*m2 + C4*m + C3
    var r =   a*m3 + C2*m2 + m + 1
  r


## Pade approximations of exp

proc expm1Pade3*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  let m2 = m*m
  splitVars:
    var xe = 0.1*m2 + 1
    var xot = (1.0/120.0)*m2 + 0.5
  let xo = m*xot
  #let n = 2*xo
  let n = xo + xo
  let d = xe - xo
  var di{.noInit.}: type(d)
  inverse(di, d)
  r := n * di
  r
proc expPade3*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  let m2 = m*m
  splitVars:
    var xe = 0.1*m2 + 1
    var xot = (1.0/120.0)*m2 + 0.5
  let xo = m*xot
  let n = xe + xo
  let d = xe - xo
  var di{.noInit.}: type(d)
  inverse(di, d)
  r := n * di
  r

# dp = v + (m'm'v-5m'vm'+vm'm')/60 + m'm'vm'm'/600
proc expm1Pade3Deriv*(m: Mat1, w: Mat2): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  let m2 = m.adj*m.adj
  splitVars:
    var xe = 0.1*m2 + 1
    var xot = (1.0/120.0)*m2 + 0.5
  let xo = m.adj*xot
  #let n = xe + xo
  let d = xe - xo
  var di{.noInit.}: type(d)
  inverse(di, d)
  let v = di*w*di
  let m2v = m2*v
  let vm2 = v*m2
  let m2vm2 = m2*vm2
  r := v
  r += (1.0/60.0)*(m2v+vm2-5*m.adj*v*m.adj)
  r += (1.0/600.0)*m2vm2
  r
proc expPade3Deriv*(m: Mat1, w: Mat2): auto {.noInit.} =
  expm1Pade3Deriv(m, w)

proc expm1Pade4*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  let m2 = m*m
  let m4 = m2*m2
  splitVars:
    var xe = (1.0/1680.0)*m4 + (3.0/28.0)*m2 + 1
    var xot = (1.0/84.0)*m2 + 0.5
  let xo = m*xot
  #let n = 2*xo
  let n = xo + xo
  let d = xe - xo
  var di{.noInit.}: type(d)
  inverse(di, d)
  r := n * di
  r
proc expPade4*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  let m2 = m*m
  let m4 = m2*m2
  splitVars:
    var xe = (1.0/1680.0)*m4 + (3.0/28.0)*m2 + 1
    var xot = (1.0/84.0)*m2 + 0.5
  let xo = m*xot
  let n = xe + xo
  let d = xe - xo
  var di{.noInit.}: type(d)
  inverse(di, d)
  r := n * di
  r

# a = m'
# v + (2aav-7ava+2vaa)/84 + a(-7aav+23ava-7vaa)a/11760 - aaavaaa/70560
proc expm1Pade4Deriv*(m: Mat1, w: Mat2): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  let m2 = m.adj*m.adj
  let m4 = m2*m2
  splitVars:
    var xe = (1.0/1680.0)*m4 + (3.0/28.0)*m2 + 1
    var xot = (1.0/84.0)*m2 + 0.5
  let xo = m.adj*xot
  #let n = xe + xo
  let d = xe - xo
  var di{.noInit.}: type(d)
  inverse(di, d)
  let v = di*w*di
  let mvm = m.adj*v*m.adj
  let m2v = m2*v
  let vm2 = v*m2
  let m2vm2 = m2*vm2
  r := v
  r += (2.0/84.0)*(m2v+vm2) + (-7.0/84.0)*mvm
  let t = (23.0/11760.0)*mvm + (-7.0/11760.0)*(m2v+vm2) + (-1.0/70560.0)*m2vm2
  r += m.adj*t*m.adj
  r
proc expPade4Deriv*(m: Mat1, w: Mat2): auto {.noInit.} =
  expm1Pade4Deriv(m, w)

proc expm1Pade5*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  let m2 = m*m
  let m4 = m2*m2
  splitVars:
    var xe = (1.0/1008.0)*m4 + (1.0/9.0)*m2 + 1
    var xot = (1.0/30240.0)*m4 + (1.0/72.0)*m2 + 0.5
  let xo = m*xot
  #let n = 2*xo
  let n = xo + xo
  let d = xe - xo
  var di{.noInit.}: type(d)
  inverse(di, d)
  r := n * di
  r
proc expPade5*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  let m2 = m*m
  let m4 = m2*m2
  splitVars:
    var xe = (1.0/1008.0)*m4 + (1.0/9.0)*m2 + 1
    var xot = (1.0/30240.0)*m4 + (1.0/72.0)*m2 + 0.5
  let xo = m*xot
  let n = xe + xo
  let d = xe - xo
  var di{.noInit.}: type(d)
  inverse(di, d)
  r := n * di
  r

proc expm1Pade6*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  let m2 = m*m
  let m4 = m2*m2
  let m6 = m2*m4
  splitVars:
    var xe = (1.0/665280.0)*m6 + (1.0/792.0)*m4 + (5.0/44.0)*m2 + 1
    var xot = (1.0/15840.0)*m4 + (1.0/66.0)*m2 + 0.5
  let xo = m*xot
  #let n = 2*xo
  let n = xo + xo
  let d = xe - xo
  var di{.noInit.}: type(d)
  inverse(di, d)
  r := n * di
  r
proc expPade6*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  let m2 = m*m
  let m4 = m2*m2
  let m6 = m2*m4
  splitVars:
    var xe = (1.0/665280.0)*m6 + (1.0/792.0)*m4 + (5.0/44.0)*m2 + 1
    var xot = (1.0/15840.0)*m4 + (1.0/66.0)*m2 + 0.5
  let xo = m*xot
  let n = xe + xo
  let d = xe - xo
  var di{.noInit.}: type(d)
  inverse(di, d)
  r := n * di
  r

proc expm1Pade7*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  let m2 = m*m
  let m4 = m2*m2
  let m6 = m2*m4
  splitVars:
    var xe = (1.0/308880.0)*m6 + (5.0/3432.0)*m4 + (3.0/26.0)*m2 + 1
    var xot = (1.0/17297280.0)*m6 + (1.0/11440.0)*m4 + (5.0/312.0)*m2 + 0.5
  let xo = m*xot
  #let n = 2*xo
  let n = xo + xo
  let d = xe - xo
  var di{.noInit.}: type(d)
  inverse(di, d)
  r := n * di
  r
proc expPade7*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  let m2 = m*m
  let m4 = m2*m2
  let m6 = m2*m4
  splitVars:
    var xe = (1.0/308880.0)*m6 + (5.0/3432.0)*m4 + (3.0/26.0)*m2 + 1
    var xot = (1.0/17297280.0)*m6 + (1.0/11440.0)*m4 + (5.0/312.0)*m2 + 0.5
  let xo = m*xot
  let n = xe + xo
  let d = xe - xo
  var di{.noInit.}: type(d)
  inverse(di, d)
  r := n * di
  r

proc expm1Pade8*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  let m2 = m*m
  let m4 = m2*m2
  let m8 = m4*m4
  let m6 = m2*m4
  splitVars:
    var xe = (1.0/518918400.0)*m8 + (1.0/205920.0)*m6 + (1.0/624.0)*m4 + (7.0/60.0)*m2 + 1
    var xot = (1.0/7207200.0)*m6 + (1.0/9360.0)*m4 + (1.0/60.0)*m2 + 0.5
  let xo = m*xot
  #let n = 2*xo
  let n = xo + xo
  let d = xe - xo
  var di{.noInit.}: type(d)
  inverse(di, d)
  r := n * di
  r
proc expPade8*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  let m2 = m*m
  let m4 = m2*m2
  let m8 = m4*m4
  let m6 = m2*m4
  splitVars:
    var xe = (1.0/518918400.0)*m8 + (1.0/205920.0)*m6 + (1.0/624.0)*m4 + (7.0/60.0)*m2 + 1
    var xot = (1.0/7207200.0)*m6 + (1.0/9360.0)*m4 + (1.0/60.0)*m2 + 0.5
  let xo = m*xot
  let n = xe + xo
  let d = xe - xo
  var di{.noInit.}: type(d)
  inverse(di, d)
  r := n * di
  r

proc expm1Pade9*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  let m2 = m*m
  let m4 = m2*m2
  let m8 = m4*m4
  let m6 = m2*m4
  splitVars:
    var xe = (1.0/196035840.0)*m8 + (1.0/159120.0)*m6 + (7.0/4080.0)*m4 + (2.0/17.0)*m2 + 1
    var xot = (1.0/17643225600.0)*m8 + (1.0/4455360.0)*m6 + (1.0/8160.0)*m4 + (7.0/408.0)*m2 + 0.5
  let xo = m*xot
  #let n = 2*xo
  let n = xo + xo
  let d = xe - xo
  var di{.noInit.}: type(d)
  inverse(di, d)
  r := n * di
  r
proc expPade9*(m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  let m2 = m*m
  let m4 = m2*m2
  let m8 = m4*m4
  let m6 = m2*m4
  splitVars:
    var xe = (1.0/196035840.0)*m8 + (1.0/159120.0)*m6 + (7.0/4080.0)*m4 + (2.0/17.0)*m2 + 1
    var xot = (1.0/17643225600.0)*m8 + (1.0/4455360.0)*m6 + (1.0/8160.0)*m4 + (7.0/408.0)*m2 + 0.5
  let xo = m*xot
  let n = xe + xo
  let d = xe - xo
  var di{.noInit.}: type(d)
  inverse(di, d)
  r := n * di
  r

proc expPade5Scale5*(m: Mat1): auto {.noInit.} =
  let ms = (1.0/32.0) * m
  let e = expPade5(ms)
  let e2 = e*e
  let e4 = e2*e2
  let e8 = e4*e4
  let e16 = e8*e8
  let e32 = e16*e16
  e32


type
  ExpKind* = enum
    ekPoly, ekPade
  ExpParam* = object
    kind*: ExpKind
    order*: int
    scale*: int
    valid*: bool
proc newExpParam*: ExpParam =
  result.kind = ekPoly
  result.order = 4
  result.scale = 20

proc expm1NoScale*[R,M:Mat1](p: var ExpParam, r: var R, m: M) =
  p.valid = true
  case p.kind
  of ekPoly:
    case p.order
    of 3: r := expm1Poly3(m)
    of 4: r := expm1Poly4(m)
    of 5: r := expm1Poly5(m)
    of 6: r := expm1Poly6(m)
    of 7: r := expm1Poly7(m)
    of 8: r := expm1Poly8(m)
    of 9: r := expm1Poly9(m)
    of 10: r := expm1Poly10(m)
    of 11: r := expm1Poly11(m)
    of 12: r := expm1Poly12(m)
    else:
      #echo "unsupported expPoly order: ", p.order
      p.valid = false
  of ekPade:
    case p.order
    of 3: r := expm1Pade3(m)
    of 4: r := expm1Pade4(m)
    of 5: r := expm1Pade5(m)
    of 6: r := expm1Pade6(m)
    of 7: r := expm1Pade7(m)
    of 8: r := expm1Pade8(m)
    of 9: r := expm1Pade9(m)
    else:
      #echo "unsupported expPade order: ", p.order
      p.valid = false

proc expNoScale*[R,M:Mat1](p: var ExpParam, r: var R, m: M) =
  p.valid = true
  case p.kind
  of ekPoly:
    case p.order
    of 3: r := expPoly3(m)
    of 4: r := expPoly4(m)
    of 5: r := expPoly5(m)
    of 6: r := expPoly6(m)
    of 7: r := expPoly7(m)
    of 8: r := expPoly8(m)
    of 9: r := expPoly9(m)
    of 10: r := expPoly10(m)
    of 11: r := expPoly11(m)
    of 12: r := expPoly12(m)
    else:
      #echo "unsupported expPoly order: ", p.order
      p.valid = false
  of ekPade:
    case p.order
    of 3: r := expPade3(m)
    of 4: r := expPade4(m)
    of 5: r := expPade5(m)
    of 6: r := expPade6(m)
    of 7: r := expPade7(m)
    of 8: r := expPade8(m)
    of 9: r := expPade9(m)
    else:
      #echo "unsupported expPade order: ", p.order
      p.valid = false

proc expm1*(p: var ExpParam, m: Mat1): auto {.noInit.} =
  var r{.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  p.valid = true
  if p.scale > 0:
    let s = numberType(m)(1.0/float(1 shl p.scale))
    var ms = s * m
    p.expm1NoScale(r, ms)
    if p.valid:
      for i in 1..p.scale:
        ms := r
        var t = ms
        t += 2
        r := ms*t
  else:
    p.expm1NoScale(r, m)
  r

proc expm1DerivNoScale*[R,M,W:Mat1](p: var ExpParam, r: var R, m: M, w: W) =
  p.valid = true
  case p.kind
  of ekPoly:
    case p.order
    of 3: r := expm1Poly3Deriv(m, w)
    of 4: r := expm1Poly4Deriv(m, w)
    #of 5: r := expm1Poly5(m)
    #of 6: r := expm1Poly6(m)
    #of 7: r := expm1Poly7(m)
    #of 8: r := expm1Poly8(m)
    #of 9: r := expm1Poly9(m)
    #of 10: r := expm1Poly10(m)
    #of 11: r := expm1Poly11(m)
    #of 12: r := expm1Poly12(m)
    else:
      #echo "unsupported expPoly order: ", p.order
      p.valid = false
  of ekPade:
    case p.order
    of 3: r := expm1Pade3Deriv(m, w)
    of 4: r := expm1Pade4Deriv(m, w)
    #of 5: r := expm1Pade5(m)
    #of 6: r := expm1Pade6(m)
    #of 7: r := expm1Pade7(m)
    #of 8: r := expm1Pade8(m)
    #of 9: r := expm1Pade9(m)
    else:
      #echo "unsupported expPade order: ", p.order
      p.valid = false

#  < 2 e' w + e' e' w >  ->  < 2 d w + d e' w + e' d w > = < d (2w+e'w+we') >
# 2(2e+ee)+(2e+ee)(2e+ee)/4 = e+3/2ee+eee+eeee/4
#  -> w + 3/2(ew+we) + (eew+ewe+wee) + (eeew+eewe+ewee+weee)/4

proc expm1Deriv*(p: var ExpParam, m: Mat1, w: Mat2): auto {.noInit.} =
  var r {.noInit.}: MatrixArray[m.nrows,m.ncols,type(m[0,0])]
  p.valid = true
  if p.scale > 0:
    let s = numberType(m)(1.0/float(1 shl p.scale))
    var ms = s * m
    var e {.noInit.}: type(r)
    p.expm1NoScale(e, ms)
    if p.valid:
      var we = 0.5*(w*e.adj + e.adj*w) + w
      for i in 2..p.scale:
        let t = e
        var t2 = e
        t2 += 2
        e := t*t2
        we += 0.5*(we*e.adj + e.adj*we)
      p.expm1DerivNoScale(r, ms, we)
  else:
    p.expm1DerivNoScale(r, m, w)
  r

proc exp*(p: var ExpParam, m: Mat1): auto {.noInit.} =
  result = p.expm1(m)
  if p.valid:
    result += 1

proc expDeriv*(p: var ExpParam, m: Mat1, w: Mat2): auto {.noInit.} =
  result = p.expm1Deriv(m, w)

when isMainModule:
  import os, strutils, strformat, fenv
  import complexNumbers
  var scale = 0
  if paramCount() > 0:
    scale = parseInt paramStr(1)
  let f32max = float32.maximumPositiveValue
  let l32max = ln(f32max)
  echo "float32.max: ", f32max, " log: ", l32max
  let f64max = float64.maximumPositiveValue
  let l64max = ln(f64max)
  echo "float64.max: ", f64max, " log: ", l64max

  proc setAH(a: var Mat1) =
    let N = a.nrows
    for i in 0..<N:
      let fi = i.float
      for j in 0..<i:
        let fj = j.float
        let tr = 0.5 + 0.7/(0.9+1.3*fi-fj)
        let ti = 0.1 + 0.3/(0.4+fi-1.1*fj)
        a[i,j].re := tr
        a[i,j].im := ti
        a[j,i].re := -tr
        a[j,i].im := ti
      let ti = 0.1 + 0.3/(0.4+fi)
      a[i,i].re := 0
      a[i,i].im := ti
    let n = 1.0/sqrt(a.norm2)
    a *= n

  proc testr() =
    var
      min = -5.0
      max = 5.0
      n = 8
      x: MatrixArray[1,1,float]
      v = newSeq[float](n)
      e = newSeq[float](n)
      p: ExpParam
      valid: bool
      nreps = 1000000
    p.scale = scale
    var s = "Kind  Ord"
    for i in 0..<n:
      let t = min + ((max-min)*i)/(n-1)
      v[i] = t
      s &= &"{t:12.4f}"
    echo s
    for kind in ExpKind:
      p.kind = kind
      for order in 3..12:
        p.order = order
        valid = false
        var secs = 0.0
        for i in 0..<n:
          x[0,0] = v[i]
          var y = p.exp(x)
          if p.valid:
            valid = true
            let ye = exp(x[0,0])
            let re = y[0,0]/ye - 1
            e[i] = re
            tic()
            for rep in 1..nreps:
              y += p.exp(x)
            secs += getElapsedTime()

        if valid:
          var o = $p.kind & " " & p.order|2
          for i in 0..<n:
            o &= &"{e[i]:12.4e}"
          echo o, "  ", ((secs*1e9)/(nreps*n))|3, " ns"

  testr()

  proc testi() =
    var
      #min = -5.0
      #max = 5.0
      min = -3.1416
      max = 3.1416
      n = 8
      x: MatrixArray[1,1,ComplexType[float]]
      v = newSeq[float](n)
      e = newSeq[float](n)
      p: ExpParam
      valid: bool
      nreps = 1000000
    p.scale = scale
    var s = "Kind  Ord"
    for i in 0..<n:
      let t = min + ((max-min)*i)/(n-1)
      v[i] = t
      s &= &"{t:12.4f}"
    echo s
    for kind in ExpKind:
      p.kind = kind
      for order in 3..12:
        p.order = order
        valid = false
        var secs = 0.0
        for i in 0..<n:
          x[0,0] := asImag(v[i])
          var y = p.exp(x)
          if p.valid:
            valid = true
            let ye = exp(x[0,0])
            #let re = sqrt((y[0,0]-ye).norm2/ye.norm2)
            let re = sqrt((y[0,0]-ye).norm2)
            e[i] = re
            tic()
            for rep in 1..nreps:
              y += p.exp(x)
            secs += getElapsedTime()

        if valid:
          var o = $p.kind & " " & p.order|2
          for i in 0..<n:
            o &= &"{e[i]:12.4e}"
          echo o, "  ", ((secs*1e9)/(nreps*n))|3, " ns"

  #testi()

  proc testm() =
    #var x: MatrixArray[3,3,float]
    var x: MatrixArray[3,3,ComplexType[float]]
    var p: ExpParam
    var nreps = 1000000
    p.scale = scale
    for kind in ExpKind:
      p.kind = kind
      for order in 3..12:
        p.order = order
        var y = p.exp(x)
        if p.valid:
          tic()
          var secs = 0.0
          for rep in 1..nreps:
            y += p.exp(x)
          secs += getElapsedTime()
          echo &"{p.kind} {p.order:-2d} {(secs*1e9)/nreps:6.1f} ns"
  #testm()

  proc testr2() =
    var
      min = -1.0
      max = 1.0
      n = 21
      x: MatrixArray[1,1,float]
    for i in 0..<n:
      x[0,0] = min + ((max-min)*i)/(n-1)
      let y3 = expPade3(x)
      let y5 = expPade5(x)
      let y7 = expPade7(x)
      let y5s = expPade5Scale5(x)
      let ye = exp(x[0,0])
      let re3 = (y3[0,0]/ye-1)
      let re5 = (y5[0,0]/ye-1)
      let re7 = (y7[0,0]/ye-1)
      let re5s = (y5s[0,0]/ye-1)
      echo x[0,0]|(-10,3), re3|(-12,4), re5|(-12,4), re7|(-12,4), re5s|(-12,4)
  #[
  proc testi() =
    var
      min = -1.0
      max = 1.0
      n = 21
      #x: MatrixArray[1,1,ImagProxy[float]]
      x: MatrixArray[1,1,ComplexType[float]]
    for i in 0..<n:
      x[0,0] = newComplex(0.0, min + ((max-min)*i)/(n-1))
      let y3 = expPade3(x)
      let y5 = expPade5(x)
      let y7 = expPade7(x)
      let y5s = expPade5Scale5(x)
      let ye = exp(x[0,0])
      let re3 = (y3[0,0]-ye)/abs(ye)
      let re5 = (y5[0,0]-ye)/abs(ye)
      let re7 = (y7[0,0]-ye)/abs(ye)
      let re5s = (y5s[0,0]-ye)/abs(ye)
      echo x[0,0]|(-12,4), re3|(-12,4), re5|(-12,4), re7|(-12,4), re5s|(-12,4)
  #testr()
  #testi()
  ]#

  proc testDeriv(T: typedesc) =
    var p: ExpParam
    var a,b,c: T
    let N = a.nrows
    let da = 1e-5
    p.scale = scale
    for kind in ExpKind:
      p.kind = kind
      for order in 3..12:
        p.order = order
        setAH(a)
        a += 0.1
        let x = p.exp(a)
        if p.valid:
          setAH(c)
          c += 1
          let f = p.expDeriv(a, c)
          if p.valid:
            echo "Testing ", N, "x", N, " ", p.kind, " ", p.order
            let tx = trace(c.adj*x).re
            for i in 0..<N:
              for j in 0..<N:
                block:
                  b := a
                  b[i,j] += da
                  var y = p.exp(b)
                  let ty = trace(c.adj*y).re
                  let d = (ty-tx)/da
                  let re = abs((d-f[i,j].re)/d)
                  if re > 4*N*da:
                    echo " ", i, " ", j, " re ", f[i,j].re, " ", d, " ", re
                block:
                  b := a
                  b[i,j] += newImag(da)
                  var y = p.exp(b)
                  let ty = trace(c.adj*y).re
                  let d = (ty-tx)/da
                  let re = abs((d-f[i,j].im)/d)
                  if re > 4*N*da:
                    echo " ", i, " ", j, " im ", f[i,j].im, " ", d, " ", re

  type
    Cmplx[T] = ComplexType[T]
    CM[N:static[int],T] = MatrixArray[N,N,Cmplx[T]]

  testDeriv(CM[1,float])
  testDeriv(CM[2,float])
  testDeriv(CM[3,float])

  # example from https://blogs.mathworks.com/cleve/2012/07/23/a-balancing-act-for-the-matrix-exponential
  proc getm3: MatrixArray[3,3,float] =
    let
      a = 2e10
      b = 4e8/6
      c = 200/3
      d = 3.0
      e = 1e-8
    result[0,1] = e
    result[1,0] = -(a+b)
    result[1,1] = -d
    result[1,2] = a
    result[2,0] = c
    result[2,2] = -c
  proc getm3exp: MatrixArray[3,3,float] =
    result[0,0] = 0.446849468283175
    result[0,1] = 1.54044157383952e-09
    result[0,2] = 0.462811453558774
    result[1,0] = -5743067.77947947
    result[1,1] = -0.0152830038686819
    result[1,2] = -4526542.71278401
    result[2,0] = 0.447722977849494
    result[2,1] = 1.54270484519591e-09
    result[2,2] = 0.463480648837651

  proc testgetm3 =
    let a = getm3()
    let ae = getm3exp()
    var p = newExpParam()
    p.scale = 60
    p.order = 12
    let b = p.exp(a)
    let d = b-ae
    echo a
    echo b
    echo sqrt(d.norm2/ae.norm2)
  testgetm3()

