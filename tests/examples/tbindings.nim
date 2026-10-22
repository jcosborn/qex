import base
import hwloc/capi
import strutils

template ifSuccess(run,cont:untyped):untyped =
  block:
    let err = run
    if err < 0:
      echo "HWLOC error: ",err," <- ",astToStr(run)
      stdout.flushFile
    else:
      cont

qexInit()
var c = getDefaultComm()
echo "rank ",myRank,"/",nRanks
threads: echo "thread ",threadNum,"/",numThreads
echo "hwloc compile-time API version: 0x", toHex(HWLOC_API_VERSION.int,8)
echo "hwloc run-time API version:     0x", toHex(hwloc_get_api_version().int,8)
const buflen = 64
var
  topology:hwloc_topology_t
  rootobj:hwloc_obj_t
  policy:hwloc_membind_policy_t

ifSuccess hwloc_topology_init(topology.addr):
  defer: hwloc_topology_destroy(topology)
  ifSuccess hwloc_topology_load(topology):
    var set = hwloc_bitmap_alloc()
    defer: hwloc_bitmap_free(set)
    #let cset_cpu = hwloc_topology_get_topology_cpuset(topology)
    let topodepth = hwloc_topology_get_depth(topology)
    echo "topology depth: ", topodepth
    let n0 = hwloc_get_nbobjs_by_depth(topology, 0)
    echo "nbobjs at depth 0: ", n0
    rootobj = hwloc_get_obj_by_depth(topology, 0, 0)
    let cset_cpu = hwloc_bitmap_dup(rootobj.cpuset);

    ifSuccess hwloc_get_cpubind(topology, set, HWLOC_CPUBIND_PROCESS.cint):
      var buffer = newStringOfCap buflen
      buffer.setlen buflen
      ifSuccess hwloc_bitmap_snprintf(cstring buffer, buflen, set):
        var
          ncur = hwloc_bitmap_weight(set)
          ntot = hwloc_bitmap_weight(cset_cpu)
        if myRank > 0:
          c.pushSend(0, buffer.cstring, buflen)
          c.pushSend(0, ncur.addr, sizeof ncur)
          c.pushSend(0, ntot.addr, sizeof ntot)
          c.waitSends
        else:
          for r in 0..<nRanks:
            if r > 0:
              c.pushRecv(r, buffer.cstring, buflen)
              c.pushRecv(r, ncur.addr, sizeof ncur)
              c.pushRecv(r, ntot.addr, sizeof ntot)
              c.waitRecvs
            echo "rank ",r," binds to cpuset ",buffer.cstring," using ",ncur," of ",ntot

    # Each thread records its binding and only the main thread communicates:
    # the comms are not thread safe, and under refc they would grow their
    # queues on a worker's heap.
    var nth = getMaxThreads()
    var tset = newSeq[array[buflen,char]](nth)
    var tcur = newSeq[cint](nth)
    threads:
      var set = hwloc_bitmap_alloc()
      defer: hwloc_bitmap_free(set)
      ifSuccess hwloc_get_cpubind(topology, set, HWLOC_CPUBIND_THREAD.cint):
        ifSuccess hwloc_bitmap_snprintf(cast[cstring](addr tset[threadNum]), buflen, set):
          tcur[threadNum] = hwloc_bitmap_weight(set)
    if myRank > 0:
      c.pushSend(0, nth.addr, sizeof nth)
      c.pushSend(0, tset[0].addr, nth*buflen)
      c.pushSend(0, tcur[0].addr, nth*sizeof(cint))
      c.waitSends
    else:
      for t in 0..<nth:
        echo "rank ",myRank," thread ",t," binds to cpuset ",cast[cstring](addr tset[t])," using ",tcur[t]
      for r in 1..<nRanks:
        c.pushRecv(r, nth.addr, sizeof nth)
        c.waitRecvs
        tset.setLen nth
        tcur.setLen nth
        c.pushRecv(r, tset[0].addr, nth*buflen)
        c.pushRecv(r, tcur[0].addr, nth*sizeof(cint))
        c.waitRecvs
        for t in 0..<nth:
          echo "rank ",r," thread ",t," binds to cpuset ",cast[cstring](addr tset[t])," using ",tcur[t]
      stdout.flushFile

  #let cset = hwloc_topology_get_topology_nodeset(topology)
  let cset = rootobj.nodeset
  var set = hwloc_bitmap_alloc()
  defer: hwloc_bitmap_free(set)
  ifSuccess hwloc_get_membind(topology, set, policy.addr, HWLOC_MEMBIND_BYNODESET.cint):
    if myRank > 0:
      c.pushSend(0, policy.addr, sizeof policy)
      c.waitSends
    else:
      for r in 0..<nRanks:
        if r > 0:
          c.pushRecv(r, policy.addr, sizeof policy)
          c.waitRecvs
        echo "rank ",r," membind policy: ",policy
    var buffer = newStringOfCap buflen
    buffer.setlen buflen
    ifSuccess hwloc_bitmap_snprintf(cstring buffer, buflen, set):
      var
        ncur = hwloc_bitmap_weight(set)
        ntot = hwloc_bitmap_weight(cset)
      if myRank > 0:
        c.pushSend(0, buffer.cstring, buflen)
        c.pushSend(0, ncur.addr, sizeof ncur)
        c.pushSend(0, ntot.addr, sizeof ntot)
        c.waitSends
      else:
        for r in 0..<nRanks:
          if r > 0:
            c.pushRecv(r, buffer.cstring, buflen)
            c.pushRecv(r, ncur.addr, sizeof ncur)
            c.pushRecv(r, ntot.addr, sizeof ntot)
            c.waitRecvs
          echo "rank ",r," binds to numa node ",buffer.cstring," using ",ncur," of ",ntot
          stdout.flushFile

qexFinalize()
