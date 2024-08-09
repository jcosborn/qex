import field, comms/comms

proc wallSource*(f: Field, t0: int, v: auto) =
  for i in f.l.sites:
    let t = f.l.coords[^1][i]
    if t == t0:
      f{i} := v

proc norm2slice*(f: SomeField, s: int): seq[float] =
  #for i in l:
  #  let t = l.vcoords[3][i][0]
  let ns = f.l.physGeom[s]
  var c = newSeq[float](ns)
  for i in f.sites:
    let k = f.l.coords[s][i]
    c[k] += f{i}.norm2
  threadRankSum c
  c
