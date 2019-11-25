import os, strformat, strutils
import qdpxxwrap, cppstring

when existsEnv("CHROMADIR"):
  const chroma_path = getEnv("CHROMADIR")
else:
  const chroma_path = getEnv("HOME")&"/lqcd/install/chroma"
const chroma_config   = chroma_path & "/bin/chroma-config"
const chroma_cxxflags = gorge(chroma_config & " --cxxflags")
const chroma_ldflags  = gorge(chroma_config & " --ldflags")
const chroma_libs     = gorge(chroma_config & " --libs")
const chroma_include  = chroma_path & "/include"

echo "\n\n\n"
echo "chroma_config= ", chroma_config
echo "cxxflags= ", chroma_cxxflags
echo "ldflags= ", chroma_ldflags
echo "libs= ", chroma_libs

{.passC: chroma_cxxflags .}
{.passL: chroma_ldflags & " " & chroma_libs .}

proc chromaInitialize*(argc: ptr cint; argv: ptr cstringArray) {.
  header: chroma_include & "/init/chroma_init.h",
  importcpp: "Chroma::initialize(@)" .}

proc chromaInitialize*() =
  echo "arg-ness"
  var argc {.importc: "cmdCount", global.}: cint
  var argv {.importc: "cmdLine", global.}: cstringArray

  echo "Initialize chroma"
  chromaInitialize(argc.addr, argv.addr)

proc chromaFinalize*() {.
  header: chroma_include & "/init/chroma_init.h",
  importcpp: "Chroma::finalize()" .}

type
  AbsInlineMeasurement* {.
    header: chroma_include & "/meas/inline/abs_inline_measurement.h",
    importcpp: "Chroma::AbsInlineMeasurement".} = object

proc chromaAbsInlineMeas*(this: AbsInlineMeasurement, update_no: int,
                          xml_out: var XMLBufferWriter) {.
  header: chroma_include & "/meas/inline/abs_inline_measurement.h",
  importcpp: "#.operator()(@)".}

proc chromaRegisterInlineAggregate*(): bool {.
  header: chroma_include & "/meas/inline/inline_aggregate.h",
  importcpp: "Chroma::InlineAggregateEnv::registerAll(@)".}

proc chromaReadInlineMeas*(xml: XmlReader, path: cstring):
  ptr AbsInlineMeasurement {.
  header: chroma_include & "/io/inline_io.h",
  importcpp: "Chroma::readInlineMeasurement(@)".}

proc chromaCreateNamedObj*(t: typedesc, s: StdString) {.
  header: "meas/inline/io/named_objmap.h",
  importcpp: "Chroma::TheNamedObjMap::Instance().create<#'1>(#)" .}

proc chromaGetNamedObj*[T](t: typedesc[T], s: StdString): var T {.
  header: "meas/inline/io/named_objmap.h",
  importcpp: "Chroma::TheNamedObjMap::Instance().getData<#'1>(#)" .}

type
  ChromaGauge = multi1d[QdpLatticeColorMatrix]

proc newChromaGauge*(): ChromaGauge =
  result.resize(1)
  #for i in 0

proc chromaCreateNamedGaugeRaw*(s: StdString) {.
  header: "meas/inline/io/named_objmap.h",
  importcpp: "Chroma::TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(#)" .}

proc chromaCreateNamedGaugeResize*(s: StdString) {.
  header: "meas/inline/io/named_objmap.h",
  importcpp: "(Chroma::TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(#)).resize(Nd)" .}

proc chromaCreateNamedGauge*(id: StdString) =
  chromaCreateNamedObj(ChromaGauge, id)
  let nd = QdpLatticeSize().size()
  chromaGetNamedObj(ChromaGauge, id).resize(nd)

proc chromaGetNamedGauge*(id: StdString): var ChromaGauge {.
  header: "meas/inline/io/named_objmap.h",
  importcpp: "Chroma::TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(#)" .}




proc test1(lat0: openarray[int]) =
  echo "Initializing Chroma ..."
  chromaInitialize()

  echo "Setting layout ..."
  QdpLayoutSetLattSize(lat0)

  echo "Creating layout ..."
  QdpLayoutCreate()

  let lat = QdpLatticeSize()
  echo lat

  var id = newStdString("test")
  chromaCreateNamedGauge(id)
  let g = addr chromaGetNamedGauge(id)
  let n = g[].size
  echo "Nd: ", n
  for i in 0..<n:
    var t = g[][i]

proc test2() =
  let id = newStdString("test")
  let g = addr chromaGetNamedGauge(id)
  let n = g[].size
  echo "Nd: ", n
  g[][0][0][0,0].real = 2.0
  g[][0][0][0,0].imag = 3.0
  let lcm = g[][0]
  let cm = lcm[0]
  let z = cm[0,0]
  let zr = z.real
  let zi = z.imag
  echo zr

proc test3() =
  let id = newStdString("test")
  let g = addr chromaGetNamedGauge(id)
  let n = g[].size
  echo "Nd: ", n
  let lcm = g[][0]
  let cm = lcm[0]
  let z = cm[0,0]
  let zr = z.real
  let zi = z.imag
  echo zr

import qex

proc test4(g: seq) =
  echo "test3"
  let lo = g[0].l
  let nd = lo.nDim
  let nc = g[0][0].nrows
  let id = newStdString("test")
  let cg = addr chromaGetNamedGauge(id)
  var coords = newSeq[cint](nd)
  for s in lo.singleSites:
    lo.coord(coords, myRank, s)
    let nn = QdpNodeNumber(coords)
    if nn != myRank:
      echo "site: ", coords
      echo "myRank: ", myRank, "  nn: ", nn
      qexAbort(-1)
    let ls = QdpGetLinearSiteIndex(coords)
    for mu in 0..<nd:
      for i in 0..<nc:
        for j in 0..<nc:
          var tr,ti: float
          assign(tr, g[mu]{s}[i,j].re)
          cg[][mu][ls][i,j].real = tr
          assign(ti, g[mu]{s}[i,j].im)
          cg[][mu][ls][i,j].imag = ti
  echo cg[][0][0][0,0].real

proc chromaPlaq(id: string): seq[float] =
  let xml = &"""
  <Plaq>
    <Name>PLAQUETTE</Name>
    <Frequency>1</Frequency>
    <Param>
      <version>2</version>
      <GaugeState>
        <Name>SIMPLE_GAUGE_STATE</Name>
        <GaugeBC>
          <Name>PERIODIC_GAUGEBC</Name>
        </GaugeBC>
      </GaugeState>
    </Param>
    <NamedObject>
      <gauge_id>{id}</gauge_id>
    </NamedObject>
  </Plaq>"""
  #echo xml
  var istr = newIStringStream(xml)
  var xmlIn = newXmlReader(istr)
  var measPtr = chromaReadInlineMeas(xmlIn, "/Plaq")
  var xmlOut: XMLBufferWriter
  chromaAbsInlineMeas(measPtr[], 0, xmlOut)
  #echo $xmlOut
  let xo = xmlOut.toXmlTree
  #echo xo
  result.newSeq(0)
  result.add parseFloat xo.child("plane_01_plaq").innerText
  result.add parseFloat xo.child("plane_02_plaq").innerText
  result.add parseFloat xo.child("plane_12_plaq").innerText
  result.add parseFloat xo.child("plane_03_plaq").innerText
  result.add parseFloat xo.child("plane_13_plaq").innerText
  result.add parseFloat xo.child("plane_23_plaq").innerText
  #echo result

when isMainModule:
  defaultSetup()
  g.random()

  test1(lat)
  test2()
  test3()

  let pl0 = plaq(g)
  echo pl0
  for i in 0..<pl0.len:
    echo i, ": ", 6.0*pl0[i]
  test4(g)

  echo chromaRegisterInlineAggregate()
  let pl1 = chromaPlaq("test")
  echo pl1

  echo "Finalizing Chroma ..."
  chromaFinalize()
  echo "Exit!"
