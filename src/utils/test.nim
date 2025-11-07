import strutils
import base

type QEXTestInfo = object
  t:float
  filename:string
  line:int
  column:int
  msg:string

type QEXTest* = ref object
  tag:string
  hidden:int
  tbegin,tend:float
  info:seq[QEXTestInfo]
  passed:bool
  subtest:seq[QEXTest]

proc newQEXTest*(tag:string, hidden=0):QEXTest =
  result.new
  result.tag = tag
  result.hidden = hidden
  result.tbegin = qexTime()
  result.tend = 0.0
  result.info.newseq(0)
  result.passed = true
  if hidden<=0:
    qexLogT result.tbegin, "QEXTest:" & tag & "."

proc newTest*(test:QEXTest, tag:string, hidden=0):QEXTest =
  result = newQEXTest(test.tag & ":" & tag, hidden+test.hidden)
  test.subtest.add(result)

proc getResult*(test:QEXTest, verbose=0):tuple[failed:int, total:int] =
  ## verbose: <hidden: fail(info), ==hidden: pass/fail(info), >hidden: always info
  ## verbose decreases by one when on subtests
  let nt = test.subtest.len
  var
    subfailed = 0
    totaltests = if nt==0: 1 else: 0    # only counts leaf subtests
    totalfailed = 0
  for i in 0..<nt:
    let (totf, tot) = test.subtest[i].getResult(verbose-1)
    if test.tend<test.subtest[i].tend:
      test.tend = test.subtest[i].tend
    totaltests += tot
    totalfailed += totf
    if totf>0:
      inc subfailed
    # echo "getResult:",test.tag," ",i,"/",nt," totf/tot ",totf,"/",tot," subfailed ",subfailed
  let
    dt = test.tend-test.tbegin
    t = if dt<1.0: formatFloat(dt*1e3, ffDecimal, 3) & " ms" else: formatFloat(dt, ffDecimal, 3) & " s"
  # echo "getResult:",test.tag," t ",test.tbegin," - ",test.tend
  if subfailed>0:
    test.passed = false
  if test.passed:
    if verbose>=test.hidden:
      if nt>0:
        echo "QEXTest: Passed. ",test.tag,":\t ",t,", ",nt," sub-tests (total ",totaltests," tests)"
      else:
        echo "QEXTest: Passed. ",test.tag,":\t ",t
  else:
    if subfailed==0:    # this is the failed leaf subtest
      inc totalfailed    # only counts subtests
      echo "QEXTest: Failed. ",test.tag,":\t ",t," ms"
  if (not test.passed) or verbose>test.hidden:
    let info = test.info
    let nlog = info.len
    if nlog>0:
      echo "    QEXTest:log:",test.tag,":"
      for j in 0..<nlog:
        echo "      [", formatFloat(info[j].t,ffDecimal,6), "s] ", info[j].filename, ":", info[j].line, ":", info[j].column, ": ", info[j].msg
  (totalfailed, totaltests)

proc qexFinalize*(test:QEXTest, verboseTest=0) =
  let (totalfailed, totaltests) = test.getResult(verboseTest)
  if totalfailed==0:
    qexFinalize()
  else:
    qexError "Failed ",totalfailed,"/",totaltests," tests"

template logInfo*(test:QEXTest, s:varargs[string,`$`]) =
  let tt = test
  let ii = instantiationInfo()
  let ti = QEXTestInfo(t:qexTime(), filename:ii.filename, line:ii.line, column:ii.column, msg:s.join)
  tt.info.add ti
  if tt.hidden<0:
    qexLogT ti.t, "QEXTest:log:", test.tag, ":", ti.filename, ":", ti.line, ":", ti.column, ": ", ti.msg

template assertAlmostEqual*(test:QEXTest, left,right:float, absTol=1e-12,relTol=1e-12) =
  let
    tt = test
    vl = left
    vr = right
    atol = absTol
    rtol = relTol
    vm = max(abs(vl), abs(vr))
    absd = abs(vl-vr)
    reld = absd/vm
    t = tt.newTest(astToStr(left)&"~"&astToStr(right), hidden=1)
  var
    res = false
    msg:string
  if absd<atol and reld<rtol:
    res = true
    msg = "Passed (absTol: " & $atol & " relTol: " & $rtol & ")"
  elif absd<atol:
    res = true    # Pass as long as absolute difference is small
    msg = "Failed relative difference (absTol: " & $atol & " relTol: " & $rtol & ")"
  elif reld<rtol:
    msg = "Failed absolute difference (absTol: " & $atol & " relTol: " & $rtol & ")"
  else:
    msg = "Failed (absTol: " & $atol & " relTol: " & $rtol & ")"
  t.logInfo vl," ~ ",vr," diff ",absd," rel ",reld," : ",msg
  t.tend = qexTime()
  t.passed = res
  tt.tend = t.tend
  tt.passed = tt.passed and res

