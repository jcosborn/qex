# Applies to thalogc.nim only: refc collects after every allocation and exports
# its object lookup (nimFulldebug), freed memory reads as 0xFF (nimBurnFree),
# and buffers handed to transfers must have a counted owner (qexBorrowCheck,
# comms/commsTypes).
switch("define", "nimFulldebug")
switch("define", "nimBurnFree")
switch("define", "qexBorrowCheck")
