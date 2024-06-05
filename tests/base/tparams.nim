import base/params

proc paramTest* =
  var i = intParam("i", 0, "Test intParam")
  var f = floatParam("f", 0.0, "Test floatParam")
  echo i, f

when isMainModule:
  import qex
  qexInit()

  var b = true
  #var x = 1

  letParam:
    bf = false
    bt = true
    bs = b
    bx = if b: true else: false
    i0 = 0
    i1 = 1
    ix = if true: 2 else: 3
    f0 = 0.0
    f1 = 1.0
    fx = if true: 2.0 else: 3.0
    s0 = "foo0"
    s1 = "foo1"
    sx = if true: "foo2" else: "foo3"
    ia0 = @[0,0,0,0]
    ia1 = @[1,1,1,1]
    iax = if true: @[2,2,2,2] else: @[3,3,3,3]
    fa0 = @[0.0,0,0,0]
    fa1 = @[1.0,1,1,1]
    fax = if true: @[2.0,2,2,2] else: @[3.0,3,3,3]

  #installLoadParams()
  #installSaveParams()
  #installHelpParam()
  installStandardParams()
  echoParams()
  processHelpParam()

  defaultSetup()
  paramTest()

  echo bf, bt, bs, bx
  echo i0, i1, ix
  echo f0, f1, fx
  echo s0, s1, sx
  echo ia0, ia1, iax
  echo fa0, fa1, fax

  processSaveParams()
  writeParamFile()
  qexFinalize()
