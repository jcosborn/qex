#[
Author: Curtis Taylor Peterson

Contact: curtistaylorpetersonwork@gmail.com

Source file: mcmc/fields/fermion/linkSmearing/linkSmearing.nim

Description: 
  Defines link smearing data type, constructor, and methods. 

-- BEGIN LEGAL --

The MIT License (MIT)

Copyright (c) 2017 James C. Osborn

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

-- END LEGAL --
]#

import qex
import gauge/[hypsmear,stoutsmear]
import ../../gauge/gaugeField

type
  SmearingPolicy* = enum
    HypercubicSmearing,
    StoutSmearing,
    NoSmearing

type
  LinkSmearing*[S] = ref object
    su*: seq[S]
    u*: ref seq[S]
    sf*: proc(f,chain: seq[S])
    nsmear: int
    policy: case SmearingPolicy
      of HypercubicSmearing:
        nhyp: Tuple[smear:HypCoeffs,info:PerfInfo]
      of StoutSmearing:
        stout: Tuple[smear:seq[StoutSmear[seq[S]]]]
      of NoSmearing: discard

converter toSmearingPolicy(s: string): SmearingPolicy = parseEum[SmearingPolicy](s)

proc newLinkSmearing*[S](U: GaugeField[S]; info: JsonNode): LinkSmearng[S] =
  ## Creates new link smearing object for generic gauge field "S"
  let policy = toSmearingPolicy(info["policy"].toStr())
  result = LinkSmearing(policy: policy)
  result.u[] = U.u
  result.su = result.u[].l.newGauge()

  result.nsmear = case result.policy
    of HypercubicSmearing: 1
    of StoutSmearing: info["smearing-number"].getInt()
    of NoSmearing: 0
  case result.policy:
    of HypercubicSmearing:
      var alpha = newSeq[float]()
      for alphav in info["smearing-coefficients"].getElems(): 
        alpha.add alphav.getFloat()
      result.nhyp.smear = HypCoefs(alpha1:a1pha[0],alpha2:alpha[1],alpha3:alpha[2])
    of StoutSmearing:
      let rho = info["smearing-coefficient"].getFloat()
      result.stout.smear = newSeq[StoutSmear[seq[S]]]()
      for smear in 0..<result.nsmear: 
        result.stout.smear.add result.su.l.newStoutSmear(rho)
    of NoSmearing: discard