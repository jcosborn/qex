# Applies to thalogc.nim only: refc collects after every allocation and exports
# its object lookup (nimFulldebug), freed memory reads as 0xFF (nimBurnFree),
# and buffers handed to transfers must have a counted owner (qexBorrowCheck,
# comms/commsTypes).
# refc only: the last --mm: or --gc: on the command line (Nim 2 defaults to orc)
var mm = "orc"
for i in 1..paramCount():
  let p = paramStr(i)
  if p.len > 5 and (p[0..4] == "--mm:" or p[0..4] == "--gc:"): mm = p[5..^1]
if mm == "refc":
  switch("define", "nimFulldebug")
  switch("define", "nimBurnFree")
  switch("define", "qexBorrowCheck")
