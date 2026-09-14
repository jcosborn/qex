import qex, gauge/plaquette

qexInit()
letParam:
  lat = latticeFromLocalLattice(@[4,4,4,4], nRanks)
  reps = 10
  trials = 5
processHelpParam()
let
  lo = lat.newLayout
  g = lo.newGauge
  f = lo.newGauge
  fp = lo.newGauge
  ds = @[lo.newGauge, lo.newGauge, lo.newGauge]
var r = lo.newRNGField(Philox4x64, 920260911'u64)
threads:
  g.random r
  for d in ds:
    d.gaussian r

template bench(label: string, body: untyped) =
  block:
    for k in 0..<2:
      body
    resetTimers()
    for trial in 0..<trials:
      tic(label)
      for _ in 0..<reps: body
      let secs = getElapsedTime()/float(reps)
      toc("sample")
      echo label, " trial=", trial, " seconds/call=", secs
    echoProf()

proc report(w: PlaqWork, label: string) =
  var nh, bytes = 0
  for row in w.h:
    for h in row:
      inc nh
      bytes += h.halo.len*sizeof(h.halo[0])
  let elem = sizeof(w.h[^1][0].halo[0]) div w.hl.lo.V
  let
    send = nh*w.hm.gather.sidx.len*elem
    recv = nh*w.hm.gather.rdest.len*elem
    local = nh*w.hm.gather.ldest.len*elem
    scratch = (w.hm.gather.sidx.len+w.hm.gather.rdest.len)*elem
  echo label, " rank=", myRank, " halo bytes=", bytes, " gathers/call=", nh,
    " send bytes/call=", send, " receive bytes/call=", recv,
    " local copy bytes/call=", local, " peak gather buffer bytes=", scratch

proc compare(label: string, x, y: auto) =
  var err: float
  threads:
    var d, n: float
    for mu in 0..<x.len:
      d += norm2(x[mu]-y[mu])
      n += norm2(x[mu])+norm2(y[mu])
    threadSingle: err = sqrt(d/(1.0+n))
  echo label, " relative difference=", err
  doAssert err < 2e-12

var total = 0.0
let ca = GaugeActionCoeffs(plaq: -float(g[0][0].nrows))
bench("production plaquette sum"):
  total += ca.gaugeAction1(g)
let wa = newPlaqWork(g[0], action=true)
report(wa, "fused plaquette sum")
bench("fused plaquette sum"):
  total += wa.plaqSum(g)
let cd = GaugeActionCoeffs(plaq: float(g[0][0].nrows))
bench("production staple sum"):
  cd.gaugeActionDeriv(g, f)
for order in 0..3:
  let
    w = newPlaqWork(g[0], order)
    seeds = ds[0..<order]
    label = "fused staple jet " & $order
  report(w, label)
  bench(label):
    w.stapleSum(g, seeds, f)
for pgm in [0.0,0.07]:
  var c = Symanzik(5.4)
  c.pgm = pgm
  let
    label = "plaq+rect pgm=" & $pgm
    neg = -1.0*c
    wp = newLoopWork(g[0])
    wf = newLoopWork(g[0])
  echo label, " coefficients: plaq=", c.plaq, " rect=", c.rect, " pgm=", c.pgm
  var ap, af: float
  bench("production " & label & " action"):
    ap = c.gaugeAction1(g,work=wp)
    total += ap
  bench("fused " & label & " action"):
    af = c.loopAction(g,work=wf)
    total += af
  let err = abs(ap-af)/(1.0+abs(ap)+abs(af))
  echo label, " action relative difference=", err
  doAssert err < 2e-12
  bench("production " & label & " force"):
    c.gaugeForce(g,fp,work=wp)
  bench("fused " & label & " force"):
    # Match gaugeForce: project U*(-grad S)^dag, including projection in timing.
    neg.loopDeriv(g,f,work=wf)
    contractProjectTAH(g,f)
  compare(label & " force",f,fp)
  forStatic order, 0, 3:
    var seeds: array[order,type(g)]
    for k in 0..<order: seeds[k] = ds[k]
    let op = if order == 0: "gradient" elif order == 1: "Hessian" else: "derivative " & $order
    when order == 1:
      bench("production " & label & " " & op):
        # gaugeDerivDeriv2 accumulates; each timed application starts from zero.
        threads:
          for mu in 0..<fp.len: fp[mu] := 0
        c.gaugeDerivDeriv2(g,ds[0],fp,work=wp)
    bench("fused " & label & " " & op):
      c.loopDeriv(g,seeds,f,work=wf)
    when order == 1: compare(label & " " & op,f,fp)
echo "checksums: ", total, " ", f.norm2
qexFinalize()
