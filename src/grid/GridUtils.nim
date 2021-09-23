import qex
import physics/qcdTypes
import grid/Grid

proc `:=`*(r: var GridLatticeGaugeField, x0: openArray[Field]) =
  var x = cast[ptr UncheckedArray[type x0[0]]](unsafeaddr x0[0])
  let lo = x[0].l
  let nd = lo.nDim
  let n = lo.nSites
  var c0 = newSeq[cint](nd)
  lo.coord(c0,lo.myrank,0)
  type Sobj = GridLatticeGaugeField.scalarObj
  var scalardata = newStdVector[Sobj](n)
  #echo lo.localGeom
  #echo c0
  threads:
    for i in lo.singleSites:
      var l = (lo.coords[^1][i].cint-c0[^1])
      for j in countdown(nd-2,0):
        l = l*lo.localGeom[j].cint + (lo.coords[j][i].cint-c0[j])
      for mu in 0..<nd:
        for ic in 0..2:
          for jc in 0..2:
            var tr,ti: float
            tr := x[mu]{i}[ic,jc].re
            ti := x[mu]{i}[ic,jc].im
            {.emit:[scalardata,"[l]._internal[mu]._internal._internal[ic][jc] = Grid::Complex(tr,ti);"].}
      #echo i, " ", l
  vectorizeFromLexOrdArray(scalardata, r)

proc `:==`*(r: var GridFermion[GridImprovedStaggeredFermionR], x: Field) =
  let lo = x.l
  let nd = lo.nDim
  let n = lo.nSites
  var c0 = newSeq[cint](nd)
  lo.coord(c0,lo.myrank,0)
  type Sobj = r.scalarObj
  var scalardata = newStdVector[Sobj](n)
  var subset = lo.getSubset("all")
  if r.checkerboard == 0:
    subset = lo.getSubset("even")
  #echo lo.localGeom
  #echo c0
  threads:
    for i in subset.singleSites:
      var l = (lo.coords[^1][i].cint-c0[^1])
      for j in countdown(nd-2,0):
        l = l*lo.localGeom[j].cint + (lo.coords[j][i].cint-c0[j])
      for ic in cint(0)..2:
        var tr,ti: float
        tr := x{i}[ic].re
        ti := x{i}[ic].im
        {.emit:[scalardata,"[l]._internal._internal._internal[ic] = Grid::Complex(tr,ti);"].}
      #echo i, " ", l
  vectorizeFromLexOrdArray(scalardata, r)

proc `:=`*(r0: var GridFermion[GridImprovedStaggeredFermionR], x: Field) =
  let r = addr r0
  let lo = x.l
  let nd = lo.nDim
  let n = lo.nSites
  var c0 = newSeq[cint](nd)
  lo.coord(c0,lo.myrank,0)
  type Sobj = r0.scalarObj
  var subset = lo.getSubset("all")
  let glsites = r0.Grid.lSites
  if glsites != n:
    if r.checkerboard == 0:
      subset = lo.getSubset("even")
    else:
      subset = lo.getSubset("odd")
  echo subset
  #threads:
  block:
    {.emit:"using namespace Grid;".}
    var t: Sobj
    #var c: Coordinate
    {.emit:"Coordinate c(4);".}
    for i in subset.singleSites:
      #echo i
      for j in 0..<nd:
        var l = lo.coords[j][i].cint - c0[j]
        {.emit:"c[j] = l;".}
      for ic in cint(0)..2:
        var tr,ti: float
        tr := x{i}[ic].re
        ti := x{i}[ic].im
        {.emit:[t,"._internal._internal._internal[ic] = Grid::Complex(tr,ti);"].}
      {.emit:["autoView(dst, ",r[],", CpuWrite);"].}
      {.emit:"pokeLocalSite(t, dst, c);".}

proc `:=`*(r0: var Field, x0: var GridFermion[GridImprovedStaggeredFermionR]) =
  let r = addr r0
  let x = addr x0
  let lo = r0.l
  let nd = lo.nDim
  let n = lo.nSites
  var c0 = newSeq[cint](nd)
  lo.coord(c0,lo.myrank,0)
  type Sobj = x0.scalarObj
  var subset = lo.getSubset("all")
  let glsites = x0.Grid.lSites
  if glsites != n:
    if x.checkerboard == 0:
      subset = lo.getSubset("even")
    else:
      subset = lo.getSubset("odd")
  echo subset
  #threads:
  block:
    {.emit:"using namespace Grid;".}
    var t: Sobj
    #var c: Coordinate
    {.emit:"Coordinate c(4);".}
    for i in subset.singleSites:
      for j in 0..<nd:
        var l = lo.coords[j][i].cint - c0[j]
        {.emit:"c[j] = l;".}
      {.emit:["autoView(dst, ",x[],", CpuRead);"].}
      {.emit:"peekLocalSite(t, dst, c);".}
      for ic in 0..2:
        var tr,ti: float
        {.emit:"tr = t._internal._internal._internal[ic].real();".}
        {.emit:"ti = t._internal._internal._internal[ic].imag();".}
        r[]{i}[ic] := newComplex(tr,ti)
