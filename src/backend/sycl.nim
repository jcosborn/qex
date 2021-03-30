import macros

{.pragma: syclh,header:"<CL/sycl.hpp>".}

type
  DefaultSelector* {.importcpp:"sycl::default_selector", syclh.} = object
  HostSelector* {.importcpp:"sycl::host_selector", syclh.} = object
  Context* {.importcpp:"sycl::context", syclh.} = object
  Device* {.importcpp:"sycl::device", syclh.} = object
  Queue* {.importcpp:"sycl::queue", syclh.} = object
  SyclBuffer*[T;N:static[int]] {.importcpp:"sycl::buffer", syclh.} = object
  AmRead {.importcpp:"sycl::access::mode::read", syclh.} = object
  AmWrite {.importcpp:"sycl::access::mode::write", syclh.} = object
  AmReadWrite {.importcpp:"sycl::access::mode::read_write", syclh.} = object
  TgHost {.importcpp:"sycl::access::target::host_buffer", syclh.} = object
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
  SyclPlus*[T] {.importcpp:"std::plus", syclh.} = object
  #SyclReduction*[T,R] {.importcpp:"sycl::intel::reduction", syclh.} = object
  #SyclRed*[N:static[int]] = object
  #SyclRed* = distinct string
  SyclRed* = distinct int

type
  Buffer*[T;N:static[int]] = ref object
    buf: SyclBuffer[T,N]
  #Accessor*[T;N:static[int];M:static[Mode];A:static[Target]] = ref object
  #  acc: SyclAccessor[T,N,M,A]
  Accessor*[T;N:static[int];M;G] = ref object
    acc: SyclAccessor[T,N,M,G]

proc newDefaultSelector*(): DefaultSelector {.importcpp:"default_selector()", syclh.}
proc newHostSelector*(): HostSelector {.importcpp:"host_selector()", syclh.}

proc selectDevice*(x: DefaultSelector): Device {.importcpp:"#.select_device()".}
proc selectDevice*(x: HostSelector): Device {.importcpp:"#.select_device()".}

type cppstring {.importcpp:"std::string",header:"string".} = object
proc `len`*(x: cppstring): cint {.importcpp:"length".}
proc cstring*(x: cppstring): cstring {.
  importcpp:"const_cast<char*>(#.c_str())".}
proc `$`*(x: cppstring): string =
  let n = x.len
  result = newString(n)
  copyMem(addr(result[0]), x.cstring, n)

proc name*(x: Device): cppstring {.
  importcpp:"#.get_info<sycl::info::device::name>()".}
proc version*(x: Device): cppstring {.
  importcpp:"#.get_info<sycl::info::device::version>()".}
proc ename*(x: Device) =
  {.emit:["printf(\"%s\\n\",", x, ".get_info<sycl::info::device::name>().c_str());"].}
proc queue*(x: Device): Queue {.importcpp:"sycl::queue(#)".}
proc maxComputeUnits*(x: Device): uint32 {.
  importcpp:"#.get_info<sycl::info::device::max_compute_units>()".}
proc preferredVectorWidthFloat*(x: Device): uint32 {.
  importcpp:"#.get_info<sycl::info::device::preferred_vector_width_float>()".}

proc wait*(q: Queue) {.importcpp:"#.wait()".}
proc device*(q: Queue): Device {.importcpp:"#.get_device()".}

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


proc mallocShared*(num_bytes: int, dev: Device, ctxt: Context):
                 pointer {.importcpp:"sycl::malloc_shared".}
proc mallocShared*[T](count: int, dev: Device, ctxt: Context):
                 ptr UncheckedArray[T] {.importcpp:"sycl::malloc_shared".}
proc mallocShared*(num_bytes: int, q: Queue):
                 pointer {.importcpp:"sycl::malloc_shared(@)".}
proc mallocShared*(T: typedesc, count: int, q: Queue): ptr UncheckedArray[T] {.
  importcpp:"sycl::malloc_shared<'1>(##,#)".}

proc syclPlus*[T](t: typedesc[T]): SyclPlus[T] {.
  importcpp:"'0()", syclh, constructor.}

#proc reduction*[T,R](x: SyclAccessor[T,1,AmReadWrite,TgGlobal], r: R):
#        SyclReduction[T,R] {.
#          importcpp:"sycl::intel::reduction(#,#)", syclh.}

#proc combine*[T,R](r: SyclReduction[T,R], x: T) {.
#  importcpp:"#.combine(#)", syclh.}

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

template submit0*(q: Queue, body: typed) =
  proc qs(qq: Queue) {.gensym.} =
    {.emit:[qq,".submit([&](cl::sycl::handler &cgh){"].}
    body
    {.emit:"});".}
  qs(q)

template submit*(q: Queue, body: typed) =
  block:
    {.emit:[q,".submit([&](cl::sycl::handler &cgh){"].}
    body
    {.emit:["});"].}

template setupSycl* =
  {.pragma: id1, importcpp:"it",nodecl,header:"",noinit,codegendecl:"".}
  {.pragma: item1, importcpp:"it",nodecl,header:"",noinit,codegendecl:"".}
setupSycl()
#macro id1*(x: untyped): untyped =
  #echo "test"
  #echo x.kind
  #echo x.repr
  #echo x.treerepr
  #x
  #newEmptyNode()

proc `[]`*(x: Id1): cint {.importcpp:"#[0]".}
proc `[]`*(x: Item1): cint {.importcpp:"#[0]".}
proc getRange*(x: Item1): cint {.importcpp:"#.get_range(0)".}

template parallelFor*(n: int, body: typed) =
  {.emit:["cgh.parallel_for<class syclkern>(sycl::range<1>{",n.uint,"},[=](sycl::item<1> it)"].}
  block:
    body
  {.emit:[");"].}

#template sum*[T](x: SyclAccessor[T,1,AmReadWrite,TgGlobal]):
#         SyclReduction[T,SyclPlus[T]] = reduction(x, syclPlus(type(T)))
var rnames {.compiletime.} = newSeq[string](0)
macro sumX*(x: typed, et: typedesc): untyped =
  let t = $x & "_red"
  let n = rnames.len
  #let et = T.getTypeInst
  #echo type(T).getImpl.repr
  #echo T.getTypeInst.repr
  #echo x.getTypeInst.repr
  echo et.repr
  rnames.add t
  #result = newNimNode(nnkPragma)
  #result.add
  result = quote do:
    {.emit:["auto ", `t`, " = sycl::ONEAPI::reduction(", `x`, ",std::plus<", `et`, ">());"].}
    #SyclRed(`t`)
    SyclRed(`n`)
    #false
  echo result.treerepr
template sum*[T](x: SyclAccessor[T,1,AmReadWrite,TgGlobal]): untyped =
  sumX(x, type(T))

macro val(x: typed): untyped =
  #echo x.repr
  let v = x.getImpl[2][1][1]
  echo v.treerepr
  result = v
macro nm(x: typed): untyped =
  let v = x.getImpl[2][1][1]
  let i = v.intVal
  echo v.treerepr
  echo i
  result = newLit(rnames[i])
macro nm2(x: typed): untyped =
  let v = x.getImpl[2][1][1]
  let i = v.intVal
  echo v.treerepr
  echo i
  result = newLit(rnames[i] & "x")
template parallelFor*(n: int, r: SyclRed, body: typed) =
  {.emit:["cgh.parallel_for(sycl::nd_range<1>{",n.uint,",1},",nm(r),",[=](sycl::nd_item<1> it,auto&",nm2(r),")"].}
  block:
    body
  {.emit:[");"].}

macro combine*(r: SyclRed, x: typed): untyped =
  echo r.repr
  echo r.getImpl.treerepr
  echo r.getImpl[2][1][1].repr
  let v = r.getImpl[2][1][1]
  let i = v.intVal
  let n = newLit(rnames[i] & "x")
  result = quote do:
    {.emit:[`n`, ".combine(", `x`, ");"].}

#template devWrite*[T;N:static[int]](b: Buffer[T,N]): untyped =
  #static: echo $T.type
  #static: echo $N
  #getAccess[T,N,amWrite,tgGlobal](b.buf)
  #b.buf.getAccess(amWrite,tgGlobal)
  #getAccess[amWrite,tgGlobal](b.buf)
  #b.buf.writeGlobal()
#  var cgh {.importcpp:"cgh",nodecl,header:"",noinit,codegendecl:"".}: Handler
#  accWrGl(b.buf, cgh)

template devRead*[T;N:static[int]](b: SyclBuffer[T,N]): untyped =
  var cgh {.importcpp:"cgh",nodecl,header:"",noinit,codegendecl:"".}: Handler
  accRdGl(b, cgh)
template devWrite*[T;N:static[int]](b: SyclBuffer[T,N]): untyped =
  var cgh {.importcpp:"cgh",nodecl,header:"",noinit,codegendecl:"".}: Handler
  accWrGl(b, cgh)
template devReadWrite*[T;N:static[int]](b: SyclBuffer[T,N]): untyped =
  var cgh {.importcpp:"cgh",nodecl,header:"",noinit,codegendecl:"".}: Handler
  accRdWrGl(b, cgh)

template submitx*(q: Queue, n: uint, body: untyped) =
  #proc job(h: Handler) = discard
  {.push stackTrace:off, lineTrace:off, line_dir:off.}
  proc qs(qq: Queue) {.gensym.} =
    {.emit:[qq,".submit([&](sycl::handler &h){"].}
    var h {.importcpp:"h",nodecl,header:"",noinit,codegendecl:"".}: Handler
    {.emit:["h.parallel_for<class axpy>(sycl::range<1>{",n,
            "}, [=] (sycl::id<1> it) {"].}
    var it {.importcpp:"it",nodecl,header:"",noinit,codegendecl:"".}: Id1
    body
    {.emit:"});});".}
  {.pop.}
  qs(q)

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

proc `[]=`*[T;G](x: SyclAccessor[T,1,AmRead,G], i: int, y: any) {.
  error:"illegal use of []= on read-only accessor".}
proc `[]=`*[T;G](x: SyclAccessor[T,1,AmRead,G], i: Id1, y: any) {.
  error:"illegal use of []= on read-only accessor".}

proc `[]=`*[T;G](x: SyclAccessor[T,1,AmWrite,G], i: int, y: any) {.
  importcpp:"#[#]=#", syclh.}
proc `[]=`*[T;G](x: SyclAccessor[T,1,AmWrite,G], i: Id1, y: any) {.
  importcpp:"#[#]=#", syclh.}

proc `[]`*[T](x: SyclHostAccessor[T,1,AmRead], i: int): T {.
  importcpp:"#[#]", syclh.}
proc `[]`*[T](x: SyclHostAccessor[T,1,AmWrite], i: int): T {.
  importcpp:"#[#]", syclh.}
proc `[]=`*[T](x: SyclHostAccessor[T,1,AmWrite], i: int, y: any) {.
  importcpp:"#[#]=#", syclh.}

when isMainModule:
  proc test =
    #let sel = DefaultSelector()
    let sel = HostSelector()
    let dev = sel.selectDevice()
    echo "Name: ", dev.name
    echo "Version: ", dev.version
    let q = dev.queue()
    let n = 100
    #let x = newBuffer[float32](n)
    let x = newSyclBuffer[float32](n)
    let y = newSyclBuffer[float32](n)
    block:
      let a = y.hostRead
      let p = a.getPointer
      for i in 0..<n:
        p[i] = i.float32
      echo p[n-1]

    q.submit:
      var xa = x.devWrite
      var ya = y.devRead
      parallelFor(n):
        var i {.id1.}: Id1
        #xa[i] = 1
        xa[i] = ya[i] + 1
    #{.pop.}
    q.wait

    block:
      let a = x.hostRead
      let p = a.getPointer
      echo p[n-1]

  test()



#[
sycl::queue q(sycl::default_selector{});

        const float A(aval);

        sycl::buffer<float,1> d_X { h_X.data(), sycl::range<1>(h_X.size()) };
        sycl::buffer<float,1> d_Y { h_Y.data(), sycl::range<1>(h_Y.size()) };
        sycl::buffer<float,1> d_Z { h_Z.data(), sycl::range<1>(h_Z.size()) };

        q.submit([&](sycl::handler& h) {

            auto X = d_X.get_access<sycl::access::mode::read>(h);
            auto Y = d_Y.get_access<sycl::access::mode::read>(h);
            auto Z = d_Z.get_access<sycl::access::mode::read_write>(h);

            h.parallel_for<class axpy>( sycl::range<1>{length}, [=] (sycl::id<1> it) {
                const int i = it[0];
                Z[i] += A * X[i] + Y[i];
            });
        });
        q.wait();

]#


