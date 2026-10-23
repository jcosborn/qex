import qex, gauge/plaquette

qexInit()
letParam:
  lat = latticeFromLocalLattice(@[4,4,4,4], nRanks)
  reps = 10
  trials = 5
  part = "all"  # hmc: only the Wilson HMC kernels
installHelpParam()
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
  # Lane gather elements are one site lane, site gather elements all lanes.
  let elem = sizeof(w.h[^1][0].halo[0]) div w.hl.lo.V
  let (gm, sm, se) = (w.hm.gather, w.hm.sites, sizeof(w.h[^1][0].halo[0]))
  let
    send = nh*(gm.sidx.len*elem+sm.sidx.len*se)
    recv = nh*(gm.rdest.len*elem+sm.rdest.len*se)
    local = nh*(gm.ldest.len*elem+sm.ldest.len*se)
    scratch = max(1, max(gm.sidx.len*elem, sm.sidx.len*se)) + max(1, max(gm.rdest.len*elem, sm.rdest.len*se))
    gathers = nh*(ord(w.hm.nsites > 0)+ord(w.hm.nlanes > 0))
  echo label, " rank=", myRank, " halo bytes=", bytes, " gathers/call=", gathers,
    " send bytes/call=", send, " receive bytes/call=", recv,
    " local copy bytes/call=", local, " required gather capacity bytes=", (if gathers > 0: scratch else: 0)

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
block:
  # The kernels of pghmc's MD step: actWilson(5.4) loop action, force and fused
  # kick, the staple force for reference, and the link update exp(a P) U with
  # and without the stored exponential. A kick at the tag of its last call keeps
  # the halo, so the halo update costs the difference of the two kick rows.
  let
    c = GaugeActionCoeffs(plaq: 5.4)
    w = newLoopWork(g[0])
    ws = newLoopWork(g[0], keep=keepStaples)
    a = -0.0125
    p = lo.newGauge
    q = lo.newGauge
    z = lo.newGauge
    e = lo.newGauge
  threads:
    p.randomTAH r
  bench("loop Wilson action"):
    total += c.loopAction(g,work=w)
  bench("loop Wilson force"):
    c.force(g,f,work=w)
  bench("staple Wilson force"):
    c.gaugeForce(g,fp,work=ws)
  compare("Wilson force",f,fp)
  bench("loop Wilson kick"):
    discard c.forceKick(g,p,q,a,work=w)
  bench("loop Wilson kick, kept halo"):
    discard c.forceKick(g,p,q,a,work=w,tag=(1,1))
  threads:
    for mu in 0..<z.len: z[mu] := a*f[mu] + p[mu]
  compare("Wilson kick",q,z)
  bench("link update"):
    threads:
      axexpahmuly(z,a,p,g)
  bench("link update with E"):
    threads:
      axexpahmuly(z,a,p,g,e)
if part == "hmc":
  echo "checksums: ", total, " ", f.norm2
  qexFinalize()
  quit()
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
for (rect, pgm) in [(false, 0.0), (true, 0.0), (true, 0.07)]:
  var c = if rect: Symanzik(5.4) else: GaugeActionCoeffs(plaq: 5.4)
  c.pgm = pgm
  let
    label = if rect: "plaq+rect pgm=" & $pgm else: "plaq"
    neg = -1.0*c
    wp = newLoopWork(g[0], keep=keepStaples)
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
