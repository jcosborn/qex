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

proc c_free(p: pointer) {.importc: "free", header: "<stdlib.h>".}

qexInit()
echo "rank ",myRank,"/",nRanks
threads: echo "thread ",threadNum,"/",numThreads
echo "hwloc compile-time API version: 0x", toHex(HWLOC_API_VERSION,8)
echo "hwloc run-time API version:     0x", toHex(hwloc_get_api_version().int,8)
var
  topology:hwloc_topology_t
  obj:hwloc_obj_t
  buffer:cstring
  policy:hwloc_membind_policy_t
ifSuccess hwloc_topology_init(topology.addr):
  defer: hwloc_topology_destroy(topology)
  ifSuccess hwloc_topology_load(topology):
    var set = hwloc_bitmap_alloc()
    defer: hwloc_bitmap_free(set)
    let cset_cpu = hwloc_topology_get_topology_cpuset(topology)
    ifSuccess hwloc_get_cpubind(topology, set, HWLOC_CPUBIND_PROCESS.cint):
      ifSuccess hwloc_bitmap_asprintf(buffer.addr, set):
        defer: c_free(buffer)
        echoAll "rank ",myrank," binds to cpu set ",hwloc_bitmap_weight(set)," / ",hwloc_bitmap_weight(cset_cpu)," PU (",buffer,")"
        stdout.flushFile

    threads:
      var buffer:cstring
      var set = hwloc_bitmap_alloc()
      defer: hwloc_bitmap_free(set)
      ifSuccess hwloc_get_cpubind(topology, set, HWLOC_CPUBIND_THREAD.cint):
        ifSuccess hwloc_bitmap_asprintf(buffer.addr, set):
          defer: c_free(buffer)
          threadCritical:
            echoAll "rank ",myrank," thread ",threadNum," binds to ",hwloc_bitmap_weight(set)," / ",hwloc_bitmap_weight(cset_cpu)," PU (",buffer,")"
            stdout.flushFile

  let cset = hwloc_topology_get_topology_nodeset(topology)
  var set = hwloc_bitmap_alloc()
  defer: hwloc_bitmap_free(set)
  ifSuccess hwloc_get_membind(topology, set, policy.addr, HWLOC_MEMBIND_BYNODESET.cint):
    echo "membind policy: ",policy
    ifSuccess hwloc_bitmap_asprintf(buffer.addr, set):
      defer: c_free(buffer)
      echoAll "rank ",myrank," binds to numa node ",hwloc_bitmap_weight(set)," / ",hwloc_bitmap_weight(cset)," numa (",buffer,")"
      stdout.flushFile
qexFinalize()
