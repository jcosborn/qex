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
      ifSuccess hwloc_bitmap_snprintf(buffer, buflen, set):
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

    threads:
      var set = hwloc_bitmap_alloc()
      defer: hwloc_bitmap_free(set)
      var buffer = newStringOfCap buflen
      buffer.setlen buflen
      ifSuccess hwloc_get_cpubind(topology, set, HWLOC_CPUBIND_THREAD.cint):
        ifSuccess hwloc_bitmap_snprintf(buffer, buflen, set):
          var
            ncur = hwloc_bitmap_weight(set)
            nth = numThreads
            tid = threadNum
          if myRank > 0:
            threadMaster:
              c.pushSend(0, nth.addr, sizeof nth)
              c.waitSends
            threadBarrier()
            threadCritical:
              c.pushSend(0, buffer.cstring, buflen)
              c.pushSend(0, ncur.addr, sizeof ncur)
              c.pushSend(0, tid.addr, sizeof tid)
              c.waitSends
          else:
            threadCritical:
              echoRaw "rank ",myRank," thread ",tid," binds to cpuset ",buffer.cstring," using ",ncur
              stdout.flushFile
            threadBarrier()
            threadMaster:
              for r in 1..<nRanks:
                c.pushRecv(r, nth.addr, sizeof nth)
                c.waitRecvs
                for t in 0..<nth:
                  c.pushRecv(r, buffer.cstring, buflen)
                  c.pushRecv(r, ncur.addr, sizeof ncur)
                  c.pushRecv(r, tid.addr, sizeof tid)
                  c.waitRecvs
                  echoRaw "rank ",r," thread ",tid," binds to cpuset ",buffer.cstring," using ",ncur
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
    ifSuccess hwloc_bitmap_snprintf(buffer, buflen, set):
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
