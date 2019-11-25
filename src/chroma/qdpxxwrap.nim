import os, cppstring
import xmlparser, xmltree
export xmltree

const qdpxxPath        = getenv("HOME")&"/lqcd/install/qdpxx"
const qdpxxInclude     = qdpxxPath & "/include"

# Multi1d

{.pragma: qdpxxMulti, header: qdpxxInclude & "/qdp_multi.h".}
type
  multi1d* {.qdpxxMulti, header: "iostream",
            importcpp: "QDP::multi1d".} [T] = object

proc constructmulti1d*[T](ns1: cint): multi1d[T] {.
  qdpxxMulti, constructor, importcpp: "QDP::multi1d<'*0>(@)".}

proc resize*[T](this: var multi1d[T]; ns1: cint) {.
  qdpxxMulti, importcpp: "#.resize(@)".}

proc size*[T](this: multi1d[T]): cint {.
  qdpxxMulti, importcpp: "#.size(@)", noSideEffect.}

proc `[]`*[T](this: multi1d[T]; i: SomeInteger): T {.
  qdpxxMulti, importcpp: "#[#]", noSideEffect.}

proc `[]=`*[T](this: var multi1d[T]; i: SomeInteger, y: T) {.
  qdpxxMulti, importcpp: "#[#] = #".}

proc `$`*(x: multi1d): string =
  let n = x.size()
  result = ""
  var sep = "["
  for i in 0..<n:
    result &= sep & $x[i]
    sep = ","
  result &= "]"


# XML support

{.pragma: qdpxml, header: qdpxxInclude & "/qdp_xmlio.h".}

type
  XmlReader* {.qdpxml, importcpp: "QDP::XMLReader".} = object
  XmlWriter* {.qdpxml, importcpp: "QDP::XMLWriter".} = object
  XmlFileWriter* {.qdpxml, importcpp: "QDP::XMLFileWriter".} = object
  XmlBufferWriter* {.qdpxml,byref,importcpp: "QDP::XMLBufferWriter".} = object

proc newXmlReader*(foo: XmlReader): XmlReader {.
  qdpxml, constructor, importcpp: "QDP::XMLReader(@)".}

proc newXmlReader*(foo: IStringStream): XmlReader {.
  qdpxml, constructor, importcpp: "QDP::XMLReader(@)".}

proc open*(this: XmlReader, foo: IStringStream) {.
  qdpxml, importcpp: "#.open(@)".}

proc str*(this: XmlBufferWriter): StdString {.
  qdpxml, importcpp: "#.str()".}

proc `$`*(this: XmlBufferWriter): string =
  $(str(this))

proc toXmlTree*(this: XmlBufferWriter): XmlNode =
  parseXml($this)

# Layout

{.pragma: qdpxx, header: qdpxxInclude & "/qdp.h".}

proc QdpLayoutSetLattSize*(nrow: multi1d[cint]) {.
  qdpxx, importcpp: "QDP::Layout::setLattSize(@)".}
  # Call the QDP layout for the given lattice size

proc QdpLayoutSetLattSize*(nrow: openArray[int]) =
  # Call the QDP layout for the given lattice size
  var jrow = constructmulti1d[cint](cint(nrow.len))
  var kk: cint = 0
  for k in items(nrow):
    jrow[kk] = cint(k)
    inc(kk)
  QDPLayoutSetLattSize(jrow)

proc QdpLayoutCreate*() {.
  qdpxx, importcpp: "QDP::Layout::create()".}

proc QdpLatticeSize*(): multi1d[cint] {.
  qdpxx, importcpp: "QDP::Layout::lattSize()".}

proc QdpGetSiteCoords*(coord: ptr cint, node,linear: cint) {.
  qdpxx, importc: "QDP::Layout::QDPXX_getSiteCoords" .}

proc QdpGetLinearSiteIndex*(coord: ptr cint): cint {.
  qdpxx, importc: "QDP::Layout::QDPXX_getLinearSiteIndex" .}

proc QdpGetLinearSiteIndex*(coord: seq[cint]): cint =
  QdpGetLinearSiteIndex(unsafeAddr coord[0])

proc QdpNodeNumber*(coord: ptr cint): cint {.
  qdpxx, importc: "QDP::Layout::QDPXX_nodeNumber" .}

proc QdpNodeNumber*(coord: seq[cint]): cint =
  QdpNodeNumber(unsafeAddr coord[0])


# Complex

type
  QdpReal* {.qdpxx, importcpp:"REAL".} = float
  QdpComplex* {.qdpxx, importcpp:"RComplex<REAL>".} = object

proc real*(x: QdpComplex): QdpReal {.
  importcpp:"#.real()".}

proc imag*(x: QdpComplex): QdpReal {.
  importcpp:"#.imag()".}

proc `real=`*(x: var QdpComplex, y: SomeNumber) {.
  importcpp:"#.real() = #".}

proc `imag=`*(x: var QdpComplex, y: SomeNumber) {.
  importcpp:"#.imag() = #".}


# Fields

type
  QdpColorMatrix* {.qdpxx, importcpp:"PColorMatrix< RComplex<REAL>, Nc>".} = object
  QdpLatticeColorMatrix* {.qdpxx, importcpp:"LatticeColorMatrix".} = object

proc `[]`*(x: QdpLatticeColorMatrix, i: SomeInteger): var QdpColorMatrix {.
  importcpp:"#.elem(#).elem()".}

proc `[]`*(x: QdpColorMatrix, i,j: SomeInteger): var QdpComplex {.
  importcpp:"#.elem(#,#)".}
