# Applies to tthreadgc.nim only: refc collects after every allocation
# (nimFulldebug) and freed memory reads as 0xFF (nimBurnFree).
import "utils.nims"
if getMm() == "refc": # refc only
  switch("define", "nimFulldebug")
  switch("define", "nimBurnFree")
