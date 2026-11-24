import macros

{.pragma: syclh,header:"<sycl/sycl.hpp>".}
{.passC: "-fsycl".}
{.passL: "-fsycl".}

type
  #DefaultSelector* {.importcpp:"sycl::default_selector_v", syclh.} = object
  #HostSelector* {.importcpp:"sycl::host_selector_v", syclh.} = object
  Context* {.importcpp:"sycl::context", syclh.} = object
  Platform* {.importcpp:"sycl::platform", syclh.} = object
  Device* {.importcpp:"sycl::device", syclh.} = object
  Queue* {.importcpp:"sycl::queue", syclh.} = object
  SyclBuffer*[T;N:static[int]] {.importcpp:"sycl::buffer", syclh.} = object
  AmRead {.importcpp:"sycl::access::mode::read", syclh.} = object
  AmWrite {.importcpp:"sycl::access::mode::write", syclh.} = object
  AmReadWrite {.importcpp:"sycl::access::mode::read_write", syclh.} = object
  #TgHost {.importcpp:"sycl::access::target::host_buffer", syclh.} = object
  TgGlobal {.importcpp:"sycl::access::target::global_buffer", syclh.} = object
  SyclAccessor*[T;N:static[int];M;A] {.
    importcpp:"sycl::accessor", syclh.} = object
  SyclHostAccessor*[T;N:static[int];M] {.
    importcpp:"sycl::host_accessor", syclh.} = object
  Handler* {.importcpp:"sycl::handler", syclh.} = object
  Id*[D:static[int]] {.importcpp:"sycl::id", syclh.} = object
  Id1* {.importcpp:"sycl::id<1>", syclh.} = object
  Item1* {.importcpp:"sycl::item<1>", syclh.} = object
  Nd1* {.importcpp:"sycl::nd_item<1>", syclh.} = object
  Group* {.importcpp:"sycl::group<1>", syclh.} = object
  SubGroup* {.importcpp:"sycl::sub_group", syclh.} = object
  WorkGroup1* {.importcpp:"sycl::work_group<1>", syclh.} = object
  LocalMultiPtr*[T] {.importcpp:"sycl::multi_ptr<'0,sycl::access::address_space::local_space,sycl::access::decorated::legacy>".} = object
  LocalPtr*[T] {.importcpp:"decltype(sycl::multi_ptr<'0,sycl::access::address_space::local_space,sycl::access::decorated::legacy>().get())",syclh.} = object
  AtomicRef*[T] {.importcpp:"sycl::atomic_ref<'0,sycl::memory_order::relaxed,sycl::memory_scope::device,sycl::access::address_space::global_space>".} = object
  SyclPlus*[T] {.importcpp:"std::plus", syclh.} = object

type
  Buffer*[T;N:static[int]] = ref object
    buf: SyclBuffer[T,N]
  #Accessor*[T;N:static[int];M:static[Mode];A:static[Target]] = ref object
  #  acc: SyclAccessor[T,N,M,A]
  Accessor*[T;N:static[int];M;G] = ref object
    acc: SyclAccessor[T,N,M,G]

#proc newDefaultSelector*(): DefaultSelector {.importcpp:"default_selector_v", syclh.}
#proc newHostSelector*(): HostSelector {.importcpp:"host_selector()", syclh.}
#var
#  defaultSelector*{.importcpp:"default_selector_v", syclh.}
#  hostSelector*{.importcpp:"host_selector_v", syclh.}

#proc selectDevice*(x: DefaultSelector): Device {.importcpp:"#.select_device()".}
#proc selectDevice*(x: HostSelector): Device {.importcpp:"#.select_device()".}
#proc selectDevice*(x: object): Device {.importcpp:"#.select_device()".}
#proc device*(x: object): Device {.importcpp:"'0(#)", constructor, syclh.}
proc defaultDevice*():Device {.importcpp:"'0{sycl::default_selector_v}", constructor, syclh.}
proc hostDevice*():Device {.importcpp:"'0{sycl::host_selector_v}", constructor, syclh.}

type cppstring {.importcpp:"std::string",header:"string".} = object
proc `len`*(x: cppstring): cint {.importcpp:"length".}
proc cstring*(x: cppstring): cstring {.
  importcpp:"const_cast<char*>(#.c_str())".}
proc `$`*(x: cppstring): string =
  let n = x.len
  result = newString(n)
  copyMem(addr(result[0]), x.cstring, n)

type cppvector[T] {.importcpp:"std::vector",header:"vector".} = object
type cppvectorIterator[T] {.importcpp:"std::vector<'0>::iterator",header:"vector".} = object
proc size*(x: cppvector): csize_t {.importcpp:"#.size()",header:"vector".}
proc `[]`*[T](x: cppvector[T], i: SomeInteger): T {.importcpp:"#[#]",header:"vector".}
proc begin*[T](x: cppvector[T]): cppvectorIterator[T] {.importcpp:"#.begin()",header:"vector".}
proc `end`*[T](x: cppvector[T]): cppvectorIterator[T] {.importcpp:"#.end()",header:"vector".}
proc max_element*[T](b,e: cppVectorIterator[T]): T {.importcpp:"*std::max_element(#,#)",header:"algorithm".}
proc max_element*(v: cppVector): auto = max_element(v.begin,v.`end`)

proc getDevices*(p: Platform): cppvector[Device] {.syclh,importcpp:"#.get_devices()".}

proc name*(x: Device): cppstring {.
  importcpp:"#.get_info<sycl::info::device::name>()".}
proc version*(x: Device): cppstring {.
  importcpp:"#.get_info<sycl::info::device::version>()".}
proc ename*(x: Device) =
  {.emit:["printf(\"%s\\n\",", x, ".get_info<sycl::info::device::name>().c_str());"].}
proc queue*(x: Device): Queue {.importcpp:"sycl::queue(#)".}
proc inOrderQueue*(x: Device): Queue {.importcpp:"sycl::queue(#, sycl::property::queue::in_order())".}
proc maxComputeUnits*(x: Device): uint32 {.
  importcpp:"#.get_info<sycl::info::device::max_compute_units>()".}
proc preferredVectorWidthFloat*(x: Device): uint32 {.
  importcpp:"#.get_info<sycl::info::device::preferred_vector_width_float>()".}
proc subGroupSizes*(x: Device): cppvector[csize_t] {.
  importcpp:"#.get_info<sycl::info::device::sub_group_sizes>()".}

proc wait*(q: Queue) {.importcpp:"#.wait()".}
proc device*(q: Queue): Device {.importcpp:"#.get_device()".}

proc getNdItem1*(): Nd1 {.importcpp:"sycl::ext::oneapi::this_work_item::get_nd_item<1>()",constructor.}
proc getWorkGroup1*(): WorkGroup1 {.importcpp:"sycl::ext::oneapi::this_work_item::get_work_group<1>()".}

proc getGlobalId1*(): csize_t {.importcpp:"sycl::ext::oneapi::this_work_item::get_nd_item<1>().get_global_id(0)".}
proc getGlobalRange1*(): csize_t {.importcpp:"sycl::ext::oneapi::this_work_item::get_nd_item<1>().get_global_range(0)".}

proc getGroup*(): Group {.importc:"sycl::ext::oneapi::this_work_item::get_work_group<1>",syclh.}
proc getSubgroup*(): SubGroup {.importc:"sycl::ext::oneapi::this_work_item::get_sub_group",syclh.}
proc size*(g: Group | SubGroup): cint {.importcpp:"#.get_local_range()[0]",syclh.}
proc range*(g: Group): cint {.importcpp:"#.get_group_range()[0]",syclh.}
proc localId*(g: Group): cint {.importcpp:"#.get_local_id()[0]",syclh.}
proc groupId*(g: Group): cint {.importcpp:"#.get_group_id()[0]",syclh.}
proc barrier*(g: Group) {.importc:"group_barrier",syclh.}
proc shift_left*[T](g: SubGroup, val: T, offset: cint):T {.importcpp:"sycl::shift_group_left(@)",syclh.}
proc localMem*[T](t: typedesc[T], g: Group):LocalMultiPtr[T] {.
  importcpp:"sycl::ext::oneapi::group_local_memory_for_overwrite<'1>(@)",syclh.}
proc get*[T](m: LocalMultiPtr[T]):LocalPtr[T] {.importcpp:"#.get()",syclh.}
proc anyOf*(g: Group, x: bool): bool {.importc:"sycl::any_of_group",syclh.}
proc makeAtomicRef*[T](x: ptr T): AtomicRef[T] {.importcpp:"'0(*#)",syclh.}
proc fetchAdd*[T](a: AtomicRef[T], val: T): T {.importcpp:"#.fetch_add(#)",syclh.}

#proc newSyclBuffer*[T](): SyclBuffer[T,0] {.noinit,
#  importcpp:"'0()", constructor, syclh.}
proc newSyclBuffer*[T](n: int): SyclBuffer[T,1] {.noinit,
  importcpp:"'0(#)", constructor, syclh.}
proc freeSyclBuffer*[T;N:static[int]](b: var SyclBuffer[T,N]) {.
  importcpp:"#.~buffer()", syclh.}
proc get_count*(b: SyclBuffer): int {.importcpp:"#.get_count()", syclh.}
template elemType*[T;N:static[int]](b: SyclBuffer[T,N]): untyped = T

#proc freeBuffer*[T;N:static[int]](x: Buffer[T,N]) =
#  freeSyclBuffer(x.buf)
#proc newBuffer*[T](n: int): Buffer[T,1] =
#  result.new(freeBuffer[T,1])
#  result.buf = newSyclBuffer[T](n)

proc newSyclAccessor*[T;N:static[int];M;A](
  b: SyclBuffer[T,N]): SyclAccessor[T,N,M,A] {.
    importcpp:"'0(#)", constructor, syclh.}
proc newSyclHostAccessor*[T;N:static[int];M](
  b: SyclBuffer[T,N]): SyclHostAccessor[T,N,M] {.
    importcpp:"'0(#)", constructor, syclh.}
proc freeSyclAccessor*[T;N:static[int];M;A](
  a: SyclAccessor[T,N,M,A]) {.importcpp:"#.~accessor()", syclh.}

proc accRdGl*[T;N:static[int]](b: SyclBuffer[T,N], cgh: Handler):
            SyclAccessor[T,N,AmRead,TgGlobal] {.
              importcpp:"'0(@)", syclh, constructor.}
proc accWrGl*[T;N:static[int]](b: SyclBuffer[T,N], cgh: Handler):
            SyclAccessor[T,N,AmWrite,TgGlobal] {.
              importcpp:"'0(@)", syclh, constructor.}
proc accRdWrGl*[T;N:static[int]](b: SyclBuffer[T,N], cgh: Handler):
              SyclAccessor[T,N,AmReadWrite,TgGlobal] {.
                importcpp:"'0(@)", syclh, constructor.}


proc freeAccessor*[T;N:static[int];M;A](
  a: Accessor[T,N,M,A]) = freeSyclAccessor(a.acc)
proc newAccessor*[T;N:static[int];M;A](
  b: Buffer[T,N]): Accessor[T,N,M,A] =
  result.new(freeAccessor[T,N,M,A])
  result.acc = newSyclAccessor[T,N,M,A](b.buf)

proc mallocDevice*(num_bytes: int, q: Queue):
    pointer {.importcpp:"sycl::malloc_device(@)".}

proc freeDevice*(p: pointer, q: Queue)
    {.importcpp:"sycl::free(@)".}

proc mallocHost*(num_bytes: int, q: Queue):
    pointer {.importcpp:"sycl::malloc_host(@)".}

proc freeHost*(p: pointer, q: Queue)
    {.importcpp:"sycl::free(@)".}

proc mallocShared*(num_bytes: int, dev: Device, ctxt: Context):
                 pointer {.importcpp:"sycl::malloc_shared".}
proc mallocShared*[T](count: int, dev: Device, ctxt: Context):
                 ptr UncheckedArray[T] {.importcpp:"sycl::malloc_shared".}
proc mallocShared*(num_bytes: int, q: Queue):
                 pointer {.importcpp:"sycl::malloc_shared(@)".}
proc mallocShared*(T: typedesc, count: int, q: Queue): ptr UncheckedArray[T] {.
  importcpp:"sycl::malloc_shared<'1>(##,#)".}

proc memcpyAsync*(q: Queue, dest,src: pointer, count: SomeInteger)
  {.importcpp:"#.memcpy(@)".}
proc memsetAsync(q: Queue, dest: pointer, val,count: SomeInteger)
  {.importcpp:"#.memset(@)".}
proc memcpy*(q: Queue, dest,src: pointer, count: SomeInteger) =
  ## waits on the queue, which covers the copy on an in-order queue
  q.memcpyAsync(dest, src, count)
  q.wait
proc memset*(q: Queue, dest: pointer, val,count: SomeInteger) =
  q.memsetAsync(dest, val, count)
  q.wait

proc syclPlus*[T](t: typedesc[T]): SyclPlus[T] {.
  importcpp:"'0()", syclh, constructor.}

#template hostRead*[T;N:static[int]](b: Buffer[T,N]): untyped =
#  newSyclAccessor[T,N,AmRead,TgHost](b.buf)
template hostRead*[T;N:static[int]](b: SyclBuffer[T,N]): untyped =
  #newSyclAccessor[T,N,AmRead,TgHost](b)
  newSyclHostAccessor[T,N,AmRead](b)
template hostWrite*[T;N:static[int]](b: SyclBuffer[T,N]): untyped =
  #newSycltAccessor[T,N,AmWrite,TgHost](b)
  newSyclHostAccessor[T,N,AmWrite](b)

#proc getPointer*[T;N:static[int];M:static[Mode];A:static[Target]](
proc getPointer*[T;N:static[int];M;A](
  a: Accessor[T,N,M,A]): ptr UncheckedArray[T] {.importcpp:"#->acc.get_pointer()", syclh.}

proc getPointer*[T;N:static[int];M;A](
  a: SyclAccessor[T,N,M,A]): ptr UncheckedArray[T] {.importcpp:"#.get_pointer()", syclh.}

proc getPointer*[T;N:static[int];M](
  a: SyclHostAccessor[T,N,M]): ptr UncheckedArray[T] {.importcpp:"#.get_pointer()", syclh.}

#template lambda(

proc `[]`*(x: Id1): cint {.importcpp:"#[0]".}
proc `[]`*(x: Item1): cint {.importcpp:"#[0]".}
proc `[]`*(x: Nd1): cint {.importcpp:"#.get_global_id(0)".}
proc getRange*(x: Item1): cint {.importcpp:"#.get_range(0)".}
proc getRange*(x: Nd1): cint {.importcpp:"#.get_global_range(0)".}

proc `[]`*[T;G](x: SyclAccessor[T,1,AmRead,G], i: int): T {.
  importcpp:"#[#]", syclh.}
proc `[]`*[T;G](x: SyclAccessor[T,1,AmRead,G], i: Id1): T {.
  importcpp:"#[#]", syclh.}
proc `[]`*[T;G](x: SyclAccessor[T,1,AmRead,G], i: Nd1): T {.
  importcpp:"#[#.get_global_id(0)]", syclh.}

proc `[]`*[T;G](x: SyclAccessor[T,1,AmWrite,G], i: int): var T {.
  importcpp:"#[#]", syclh.}
proc `[]`*[T;G](x: SyclAccessor[T,1,AmWrite,G], i: Id1): var T {.
  importcpp:"#[#]", syclh.}

proc `[]=`*[T;G](x: SyclAccessor[T,1,AmRead,G], i: int, y: auto) {.
  error:"illegal use of []= on read-only accessor".}
proc `[]=`*[T;G](x: SyclAccessor[T,1,AmRead,G], i: Id1, y: auto) {.
  error:"illegal use of []= on read-only accessor".}

proc `[]=`*[T;G](x: SyclAccessor[T,1,AmWrite,G], i: int, y: auto) {.
  importcpp:"#[#]=#", syclh.}
proc `[]=`*[T;G](x: SyclAccessor[T,1,AmWrite,G], i: Id1, y: auto) {.
  importcpp:"#[#]=#", syclh.}

proc `[]`*[T](x: SyclHostAccessor[T,1,AmRead], i: int): T {.
  importcpp:"#[#]", syclh.}
proc `[]`*[T](x: SyclHostAccessor[T,1,AmWrite], i: int): T {.
  importcpp:"#[#]", syclh.}
proc `[]=`*[T](x: SyclHostAccessor[T,1,AmWrite], i: int, y: auto) {.
  importcpp:"#[#]=#", syclh.}
