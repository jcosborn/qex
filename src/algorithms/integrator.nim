import base
import mdevolve
import strutils

export mdevolve

type IntegratorProc* = proc(T,V:Integrator; steps:int):Integrator
converter toIntegratorProc*(s:string):IntegratorProc =
  template mkProc1(s:untyped):IntegratorProc =
    proc mkInt(T,V:Integrator; steps:int):Integrator {.gensym.} =
      `mk s`(T = T, V = V, steps = steps)
    mkInt
  template mkProc2(s:untyped):IntegratorProc =
    proc mkInt(T,V:Integrator; steps:int):Integrator {.gensym.} =
      `mk s`(T = T, V = V, steps = steps, ss[1].parseFloat)
    mkInt
  template mkProc3(s:untyped):IntegratorProc =
    proc mkInt(T,V:Integrator; steps:int):Integrator {.gensym.} =
      `mk s`(T = T, V = V, steps = steps, ss[1].parseFloat, ss[2].parseFloat)
    mkInt
  template mkProc4(s:untyped):IntegratorProc =
    proc mkInt(T,V:Integrator; steps:int):Integrator {.gensym.} =
      `mk s`(T = T, V = V, steps = steps, ss[1].parseFloat, ss[2].parseFloat, ss[3].parseFloat)
    mkInt
  template mkProc5(s:untyped):IntegratorProc =
    proc mkInt(T,V:Integrator; steps:int):Integrator {.gensym.} =
      `mk s`(T = T, V = V, steps = steps, ss[1].parseFloat, ss[2].parseFloat, ss[3].parseFloat, ss[4].parseFloat)
    mkInt
  let ss = s.split(',')
  # Omelyan's triple star integrators, see Omelyan et. al. (2003)
  case ss[0]:
  of "2MN":
    if ss.len == 1: return mkProc1(Omelyan2MN)
    else: return mkProc2(Omelyan2MN)
  of "4MN5FP":
    if ss.len == 1: return mkProc1(Omelyan4MN5FP)
    elif ss.len == 2: return mkProc2(Omelyan4MN5FP)
    elif ss.len == 3: return mkProc3(Omelyan4MN5FP)
    elif ss.len == 4: return mkProc4(Omelyan4MN5FP)
    elif ss.len == 5: return mkProc5(Omelyan4MN5FP)
    else: return mkProc2(Omelyan4MN5FP)
  of "4MN5FV":
    if ss.len == 1: return mkProc1(Omelyan4MN5FV)
    elif ss.len == 2: return mkProc2(Omelyan4MN5FV)
    elif ss.len == 3: return mkProc3(Omelyan4MN5FV)
    elif ss.len == 4: return mkProc4(Omelyan4MN5FV)
    elif ss.len == 5: return mkProc5(Omelyan4MN5FV)
    else: return mkProc2(Omelyan4MN5FV)
  of "6MN7FV": return mkProc1(Omelyan6MN7FV)
  of "4MN3F1GP":  # lambda = 0.2725431326761773  is  FUEL f3g a0=0.109
    if ss.len == 1: return mkProc1(Omelyan4MN3F1GP)
    else: return mkProc2(Omelyan4MN3F1GP)
  of "4MN4F2GVG": return mkProc1(Omelyan4MN4F2GVG)
  of "4MN4F2GV": return mkProc1(Omelyan4MN4F2GV)
  of "4MN5F1GV": return mkProc1(Omelyan4MN5F1GV)
  of "4MN5F1GP": return mkProc1(Omelyan4MN5F1GP)
  of "4MN5F2GV": return mkProc1(Omelyan4MN5F2GV)
  of "4MN5F2GP": return mkProc1(Omelyan4MN5F2GP)
  of "6MN5F3GP": return mkProc1(Omelyan6MN5F3GP)
  else:
    qexError "Cannot parse integrator: '", s, "'\n",
      """Available integrators (with default parameters):
      2MN,0.1931833275037836
      4MN5FP,0.2750081212332419,-0.1347950099106792,-0.08442961950707149,0.3549000571574260
      4MN5FV,0.2539785108410595,-0.03230286765269967,0.08398315262876693,0.6822365335719091
      6MN7FV
      4MN3F1GP,0.2470939580390842
      4MN4F2GVG
      4MN4F2GV
      4MN5F1GV
      4MN5F1GP
      4MN5F2GV
      4MN5F2GP
      6MN5F3GP"""
