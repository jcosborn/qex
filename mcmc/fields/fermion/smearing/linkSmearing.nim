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

import smearingDefinitions
import nhypSmearing
import stoutSmearing
import json

type
  LinkSmearing*[N:static[int]] = ref object
    ## LinkSmearing object
    ## Attributes:
    ##   nsmear: Number of smearings; applies to only stout for now
    ##   su: Smeared gauge field (../../latticeFields.nim)
    ##   u: Reference to unsmeared gauge field (../../latticeFields.nim)
    ##   derivative: Procedure implementing chain rule
    ##      Inputs: 
    ##        f: Output of function of links
    ##        chain: df/dV for V the smeared link 
    ##    smearingPolicy: SmearingPolicy specifying which
    ##      smearing is to be implemented
    ##    nhyp: Tuple containing HypCoefs nHYP container
    ##      and container of nHYP runtime information
    ##      (src/gauge/hypsmear.nim)
    ##    stout: Sequence of "Stout" objects
    nsmear: int
    su*: seq[DComplexMatrixV[N]]
    u*: ref seq[DComplexMatrixV[N]]
    derivative*: proc(f,chain:seq[DComplexMatrixV[N]])
    case smearingPolicy: SmearingPolicy
      of HypercubicSmearing:
        nhyp: tuple[parameters:HypCoefs,info:PerfInfo]
      of StoutSmearing:
        stout: seq[Stout[N]]
      of NoSmearing: discard

proc newLinkSmearing*[N:static[int]](
    u: seq[DComplexMatrixV[N]]; 
    info: JsonNode
  ): LinkSmearing[N] =
  ## LinkSmearing object constructor
  ## Inputs:
  ##   u: Gauge field data type (../../latticeFields.nim)
  ##   info: JSON data type containing gauge field info
  ## Output:
  ##   result: Link smearing object 
  
  # Set policy variable and pick smearing from available options
  var policy: SmearingPolicy
  if not info.hasKey("smearing"): qexError SmearingError1
  case info["smearing"].getStr():
    of "nhyp","nHYP","hypercubic","Hypercubic": policy = HypercubicSmearing
    of "stout","Stout": policy = StoutSmearing
    of "none","NONE","na","NA","n/a","N/A": policy = NoSmearing
    else:
      echo SmearingError2
      throwError info["smearing"].getStr() & " is not a valid smearing option"

  # Construct link smearing object
  result = LinkSmearing[N](smearingPolicy: policy)
  
  # Create reference to unsmeared gauge field & construct smeared
  # gauge field
  new(result.u)
  result.u[] = u
  result.su = u[0].l.newComplexGaugeLinks(N)

  # Construct smearing according to chosen smearing policy
  result.nsmear = case result.smearingPolicy
    of HypercubicSmearing: 1
    of StoutSmearing: 
      (if info.hasKey("smearing-levels"): info["smearing-levels"].getInt() else: 1)
    of NoSmearing: 0
  case result.smearingPolicy:
    of HypercubicSmearing:
      var alph: seq[float]
      case info.hasKey("smearing-coefficients"):
        of true:
          alph = newSeq[float]()
          for alphav in info["smearing-coefficients"].getElems(): 
            alph.add alphav.getFloat()
        of false: alph = @[0.4,0.5,0.5] # Default if not provided
      result.nhyp.parameters = HypCoefs(alpha1:alph[0],alpha2:alph[1],alpha3:alph[2])
      qexLog result.smearingPolicy.instantiationLogMessage(alph,result.nsmear)
    of StoutSmearing:
      let rho = case info.hasKey("smearing-coefficient")
        of true: info["smearing-coefficient"].getFloat()
        of false: 0.1 # Default if not proided
      var rhos: seq[float]
      result.stout = newSeq[Stout[N]]()
      case info.hasKey("smearing-coefficients"):
        of true:
          for elem in info["smearing-coefficients"].getElems():
            rhos.add elem.getFloat()
        of false: # Default to same rho for all smears
          for smear in 0..<result.nsmear: rhos.add rho
      if not (rhos.len == result.nsmear):
        throwError "\"# stout\" != \"# rho values provided\""
      for smear in 0..<result.nsmear: result.stout.add u.newStout(rhos[smear])
      qexLog result.smearingPolicy.instantiationLogMessage(rhos,result.nsmear)
    of NoSmearing: qexLog result.smearingPolicy.instantiationLogMessage(@[],0)

proc smear*(self: LinkSmearing) =
  ## Smears gauge links
  ## Inputs:
  ##   self: LinkSmearing object
  case self.smearingPolicy:
    of HypercubicSmearing: self.nhyp.smear(self.u[],self.su)
    of StoutSmearing: self.stout.smear(self.u[],self.su)
    of NoSmearing: self.su := self.u[]

proc smearGetDerivative*(self: LinkSmearing) =
  ## Smears gauge links & returns proc for derivative of function of smeared links
  ## Inputs:
  ##   self: LinkSmearing object
  self.derivative = case self.smearingPolicy
    of HypercubicSmearing: self.nhyp.smearGetDerivative(self.u[],self.su)
    of StoutSmearing: self.stout.smearGetDerivative(self.u[],self.su)
    of NoSmearing: identity

if isMainModule:
  qexInit()

  # Create lattice
  var
    lat = @[8,8,8,8]
    lo = newLayout(lat)
    testGC = GaugeActionCoeffs(plaq: 1.0)

  # Create smearing objects
  var 
    u = lo.newComplexGaugeLinks(3)
    f = lo.newComplexGaugeLinks(3)
    nhypInfo = %* {
      "smearing": "nhyp",
      "smearing-coefficients": [0.4,0.5,0.5]
    }
    stoutInfo = %* {
      "smearing": "stout",
      "smearing-levels": 3,
      "smearing-coefficient": 0.1, # Won't be used because vvvv
      "smearing-coefficients": [0.1,0.2,0.3] # Set custom values for "rho"
    }
    noneInfo = %* {"smearing": "none"}
    nhyp = u.newLinkSmearing(nhypInfo)
    stout = u.newLinkSmearing(stoutInfo)
    none = u.newLinkSmearing(noneInfo)

  # Set gauge to random; should update reference in smearing objects
  u.random

  # Test out smearing
  nhyp.smear()
  stout.smear()
  none.smear()

  # Print out actions w/ unsmeared/smeared links
  echo "nHYP: ", testGC.gaugeAction1(nhyp.u[]), "/", testGC.gaugeAction1(nhyp.su)
  echo "stout: ", testGC.gaugeAction1(stout.u[]), "/", testGC.gaugeAction1(stout.su)
  echo "none: ", testGC.gaugeAction1(none.u[]), "/", testGC.gaugeAction1(none.su)

  # Test out smear with objective of getting force
  nhyp.smearGetDerivative()
  stout.smearGetDerivative()
  none.smearGetDerivative()

  # Smear force
  proc testForceSmear(self: LinkSmearing) =
    var f2u,f2s: float
    testGC.gaugeActionDeriv(self.u[],f)
    f2u = f.norm2
    testGC.gaugeActionDeriv(self.su,f)
    f2s = f.norm2
    self.derivative(f,f)
    echo $(self.smearingPolicy) & " smeared force f^2: ", f.norm2, "/", f2s, "/", f2u
  nhyp.testForceSmear()
  stout.testForceSmear()
  none.testForceSmear()

  qexFinalize()