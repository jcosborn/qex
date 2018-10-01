import qex

proc f2cMap(c: var seq[int32]; f: seq[int32]; lsc,lsf: seq[int]) =
  for i in 0..<c.len:
    c[i] = (f[i]*lsc[i].int32) div lsf[i].int32

proc c2fMap(fmin,fmax: var seq[int32]; c: seq[int32]; lsf,lsc: seq[int]) =
  let nd = c.len
  for i in 0..<nd:
    fmin[i] = (c[i]*lsf[i].int32+lsc[i].int32-1) div lsc[i].int32
    fmax[i] = (c[i]*lsf[i].int32+lsf[i].int32-1) div lsc[i].int32
  var c2 = newSeq[int32](nd)
  var f2 = newSeq[int32](nd)
  template CHKEQ =
    for k in 0..<nd:
      if c[k]!=c2[k]:
        echo "error c2fMap: ", k, " ", c[k].int, " ", c2[k].int
        echo "c: ", @c
        echo "fmin: ", @fmin
        echo "fmax: ", @fmax
        quit -1
  f2cMap(c2, fmin, lsc, lsf)
  CHKEQ
  f2cMap(c2, fmax, lsc, lsf)
  CHKEQ

proc checkLocal(f: Layout, c: Layout2) =
  let nd = f.nDim
  let nsc = c.nSites
  let lsf = f.physGeom
  let lsc = c.physGeom
  var xc = newSeq[int32](nd)
  var xf0 = newSeq[int32](nd)
  var xf1 = newSeq[int32](nd)
  var nonlocal = 0.0
  for i in 0..<nsc:
    for j in 0..<nd:
      xc[j] = c.coords[j][i].int32
    c2fMap(xf0, xf1, xc, lsf, lsc)
    if f.rankIndex(xf0).rank != myRank or
       f.rankIndex(xf1).rank != myRank:
      nonlocal = 1.0
      break
  rankSum(nonlocal)
  if nonLocal > 0.0:
    echo nonlocal.int
    quit -1

proc getpar[T](x: openArray[T]): T =
  for i in 0..<x.len:
    result += x[i]
  result mod 2

type
  MgBlockSites* = object
    sites: ptr carray[int32]
    nsites: int
  MgBlock*[VF,VC: static[int]] = object
    fine*: Layout[VF]
    coarse*: Layout[VC]
    lb: seq[array[2,MgBlockSites]]
    sites: seq[int32]
    csites*: seq[int32]

proc newMgBlock*[VF,VC](f: Layout[VF], c: Layout[VC]): MgBlock[VF,VC] =
  result.fine = f
  result.coarse = c
  let nd = f.nDim
  let lsf = f.physGeom
  let lsc = c.physGeom

  checkLocal(f, c)

  let nlb = c.nSites
  result.lb.newSeq(nlb)
  let nsf = f.nSites
  result.sites.newSeq(nsf)
  var xf = newSeq[int32](nd)
  var xc = newSeq[int32](nd)
  result.csites.newSeq(nsf)

  for j in 0..<nsf:
    f.coord(xf, (myRank,j))
    f2cMap(xc, xf, lsc, lsf)
    let b = c.rankIndex(xc).index
    let par = getpar(xf)
    inc result.lb[b][par].nsites
    result.csites[j] = b.int32

  var n = 0
  for b in 0..<nlb:
    for par in 0..1:
      result.lb[b][par].sites = cast[ptr cArray[int32]](addr result.sites[n])
      n += result.lb[b][par].nsites
      result.lb[b][par].nsites = 0

  for j in 0..<nsf:
    f.coord(xf, (myRank,j))
    f2cMap(xc, xf, lsc, lsf)
    let b = c.rankIndex(xc).index
    let par = getpar(xf)
    let s = result.lb[b][par].nsites
    result.lb[b][par].sites[s] = j.int32
    result.lb[b][par].nsites = s+1

when isMainModule:
  qexInit()
  let latF = [16,16,16,16]
  var loF = newLayout(latF)

  let latC = [8,8,8,8]
  var loC = newLayout(latC, loF.V, loF.rankGeom, loF.innerGeom)

  let b = newMgBlock(loF, loC)
  qexFinalize()
