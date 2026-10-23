import qex
import base/[alignedMem, profile]
import gauge/staples
import gauge/plaquette
import comms/halo
import std/[strutils, sequtils]

qexInit()
letParam:
  loc = @[4,4,4,4]  # local lattice
  lat = latticeFromLocalLattice(loc, nRanks)
  reps = 512
  trials = 5
  keep = "staples"  # LoopWork retention: none, staples
  families = "Wilson,Symanzik,Iwasaki,DBW2,adjoint"
  timers = 0  # print the timer tree of each family and operation
  path = "staple"  # staple: gaugeAction1/gaugeForce; loop: loopAction/loopDeriv; plaq: plaqSum/stapleSum (Wilson)
  reuse = 0  # staple and loop paths: the pair operation tags the links so the force reuses the action's staples or halo
installHelpParam()
processHelpParam()
if path notin ["staple", "loop", "plaq"]:
  raise newException(ValueError, "path must be staple, loop, or plaq")
let
  lo = lat.newLayout
  cases = [
    ("Wilson", GaugeActionCoeffs(plaq: 5.4)),
    ("Symanzik", Symanzik(5.4)),
    ("Iwasaki", Iwasaki(5.4)),
    ("DBW2", DBW2(5.4)),
    ("adjoint", GaugeActionCoeffs(plaq: 5.4, adjplaq: 0.2))]
var
  rng = lo.newRNGField(Philox4x64, 129911u64)
  g = lo.newGauge
  f = lo.newGauge
  total = 0.0
let work = newLoopWork(g[0], parseEnum[LoopKeep](keep))
let nd = g.len
let nc = g[0][0].ncols
threads:
  g.random rng
setRawMemGcThreshold(high(int))

proc shiftRow(len: int): tuple[send, recv, pack, msgs: int] =
  ## Buffer bytes and messages of one staple-path shift row along +-nu.
  for row in newLinkShifts(lo, typeof(g[0][0]), len):
    for s in row:
      if s.sb != nil:
        result.send += s.sb.sq.sbufSize
        result.recv += s.sb.sq.rbufSize
        result.pack += s.sb.lbufSize
        result.msgs += s.si.nSendRanks + s.si.nRecvRanks

proc haloReport(label: string, maps: openArray[HaloMap], kept: int, elem: int) =
  ## One update of a field through each map. Lane gather elements are one site
  ## lane, site gather elements all lanes.
  var send, recv, local, msgs = 0
  for hm in maps:
    let (gm, sm, se) = (hm.gather, hm.sites, elem*lo.V)
    send += gm.sidx.len*elem+sm.sidx.len*se
    recv += gm.rdest.len*elem+sm.rdest.len*se
    local += gm.ldest.len*elem+sm.ldest.len*se
    msgs += gm.smsginfo.len+gm.rmsginfo.len+sm.smsginfo.len+sm.rmsginfo.len
  echo label, " halo kept bytes ", kept, " = ", kept.float/g[0].s.mem.size.float, " fields; per call: updates ", maps.len,
    " send ", send, " recv ", recv, " local copy ", local, " bytes, messages ", msgs

proc loopReport(w: LoopWork, label: string) =
  ## Halo rows of a loop workspace: kept bytes and the exchange of one update of every slot-0 field.
  for row in w.halos:
    var kept = 0
    for slot in row.fields:
      for h in slot:
        if h != nil:
          kept += h.halo.len*sizeof(h.halo[0])
    haloReport(label, row.maps, kept, sizeof(row.fields[0][0].halo[0]) div lo.V)

let fieldBytes = g[0].s.mem.size
let fwd = shiftRow(1)
let bck = shiftRow(-1)
echo "ranks ", nRanks, " local ", lo.localGeom, " field bytes ", fieldBytes, " keep ", keep, " path ", path, " reuse ", reuse
echo "shift row +nu: send ", fwd.send, " recv ", fwd.recv, " pack ", fwd.pack, " bytes, messages ", fwd.msgs,
  "; -nu: send ", bck.send, " recv ", bck.recv, " pack ", bck.pack, " bytes, messages ", bck.msgs
echo "staple path exchange bytes/call (send+recv): action plaq ", fwd.send+fwd.recv, " rect ", fwd.send+fwd.recv+bck.send+bck.recv,
  "; force plaq ", fwd.send+fwd.recv+bck.send+bck.recv, " rect ", 2*(fwd.send+fwd.recv+bck.send+bck.recv)
let wl = newLoopWork(g[0])  # loop path: halos persist, no staples
var wa, wd: typeof(newPlaqWork(g[0]))
if path == "plaq":
  wa = newPlaqWork(g[0], action=true)
  wd = newPlaqWork(g[0], 0)
  let elem = sizeof(wa.h[0][0].halo[0]) div lo.V
  haloReport("plaq action", newSeqWith(nd, wa.hm), nd*wa.h[0][0].halo.len*sizeof(wa.h[0][0].halo[0]), elem)
  haloReport("plaq force", newSeqWith(nd, wd.hm), nd*wd.h[0][0].halo.len*sizeof(wd.h[0][0].halo[0]), elem)
echo "family,operation,trial,seconds_per_call,allocated_bytes_per_call"
var step = 0
for (name,c) in cases:
  if name notin families.split(','): continue
  if path == "plaq" and name != "Wilson": continue
  let neg = -1.0*c
  let ds0: seq[seq[typeof(g[0])]] = @[]
  var act: float
  template call(isForce: bool, t: LoopTag) =
    if c.adjplaq != 0:
      if isForce: c.forceA(g, f)
      else: act = c.actionA(g)
    elif path == "staple":
      if isForce: c.gaugeForce(g, f, work=work, tag=t)
      else: act = c.gaugeAction1(g, work=work, tag=t)
    elif path == "loop":
      if isForce:
        neg.loopDeriv(g, f, work=wl, tag=t)
        contractProjectTAH(g, f)
      else: act = c.loopAction(g, work=wl, tag=t)
    else:
      if isForce:
        wd.stapleSum(g, ds0, f)
        threads:
          for mu in 0..<nd: f[mu] *= c.plaq/float(nc)
        contractProjectTAH(g, f)
      else: act = -c.plaq/float(nc)*wa.plaqSum(g)
    if not isForce: total += act
  block:
    let want = lo.newGauge
    let expected = if c.adjplaq != 0: c.actionA(g) else: c.gaugeAction2(g)
    if c.adjplaq != 0: c.forceA(g, want)
    else: c.gaugeForce(g, want)
    call(false, (0, 0))
    let ae = abs(act-expected)/(1.0+abs(act)+abs(expected))
    call(true, (0, 0))
    var fe: float
    threads:
      var d, n: float
      for mu in 0..<nd:
        d += norm2(f[mu]-want[mu])
        n += norm2(f[mu])+norm2(want[mu])
      threadSingle: fe = sqrt(d/(1.0+n))
    echo name, " relative differences: action ", ae, " force ", fe
    doAssert ae < 2e-12 and fe < 2e-12
    for x in want: x.release
  for op in ["action", "force", "pair"]:
    if op == "pair" and path == "plaq": continue
    template step =
      if op == "pair":
        inc step
        let tag: LoopTag = if reuse != 0: (1, step) else: (0, 0)
        call(false, tag)
        call(true, tag)
      else: call(op == "force", (0, 0))
    for _ in 0..<4: step
    resetTimers()
    for trial in 0..<trials:
      GC_fullCollect()
      let before = getRawMemAllocated()
      tic("gauge action sample")
      for _ in 0..<reps: step
      let secs = getElapsedTime()/float(reps)
      toc("sample")
      let bytes = float(getRawMemAllocated()-before)/float(reps)
      echo name, ",", op, ",", trial, ",", secs, ",", bytes
    let w = if path == "loop": wl else: work
    var owned = w.bytes
    if path == "plaq":
      owned = 0
      for p in [wa, wd]:
        for row in p.h:
          for h in row: owned += h.halo.mem.size
    let buf = getDefaultComm().gatherBuffers
    let capacity = if buf == nil: 0 else: buf.send.n+buf.recv.n
    echo name, " ", op, " workspace bytes ", owned, " = ", owned.float/fieldBytes.float,
      " fields; gather capacity bytes ", capacity, "; peak RSS MiB ", maxRssMiB()
    if path != "plaq":
      echo name, " ", op, " builds ", w.builds, " reuses ", w.reuses, " halo updates ", w.updates
    if path == "loop": loopReport(wl, name & " " & op)
    if timers != 0:
      echo "timers ", name, " ", op, " calls ", trials*reps
      echoTimers()
echo "action checksum: ", total
qexFinalize()
