# utility routines for nimscript configuration files

proc getMm*(): string =
  # return GC type as string
  result = "refc"  # current QEX default is refc
  # find the last --mm: or --gc: on the command line
  for i in 1..paramCount():
    let p = paramStr(i)
    if p.len > 5 and (p[0..4] == "--mm:" or p[0..4] == "--gc:"): result = p[5..^1]
