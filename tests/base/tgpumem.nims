# Applies to tgpumem.nim only: the CPU backend keeps device copies apart from
# host memory (beSharedMem=false), refc collects after every allocation
# (nimFulldebug), and freed memory reads as 0xFF (nimBurnFree).
import "utils.nims"
switch("define", "beSharedMem=false")
if getMm() == "refc": # refc only
  switch("define", "nimFulldebug")
  switch("define", "nimBurnFree")
