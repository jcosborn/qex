#[
Author: Curtis Taylor Peterson

Contact: curtistaylorpetersonwork@gmail.com

Source file: mcmc/fields/gauge/unitary/unitary.nim

Description: 
  Defines unitary gauge field objects, constructors, and methods

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

import ../gaugeFields

type
  SpecialUnitaryGauge*[N:static[int]] = ref object of LatticeField
    ## SpecialUnitaryGauge object:
    ## Attributes:
    ##   u: gauge field (mcmc/fields/latticeFields.nim)
    ##   actionPolicy: named gauge action (mcmc/fields/gauge/gaugeFields.nim)
    ##   actionParams: gauge action coefficients (src/gauge/gaugeAction.nim)
    u*: seq[DComplexMatrixV[N:static[int]]] 
    actionPolicy: GaugeActionPolicy 
    actionParams: GaugeActionCoeffs 

# Main "SpecialUnitaryGauge" constructor
proc newSpecialUnitaryGauge(
    l: Layout; 
    info: JsonNode;
    n: static[int];
  ): auto =
  ## Creates new "SpecialUnitaryGauge" object
  ## Inputs:
  ##   l: lattice layout (src/layout/)
  ##   info: JsonNode (species information about U1 or SU(N) gauge field)
  ##   n: number of "colors"
  ## Output:
  ##   result: SpecialUnitaryGauge object
  
  # Throw error if action or bare gauge coupling not specified; 
  # otherwise, save temporarily
  if not info.hasKey("action"): qexError GaugeError1
  if not info.hasKey("beta"): qexError GaugeError2
  let
    actionPolicy = toGaugeGroupPolicy(info["action"].getStr())
    beta = info["beta"].getFloat()

  # Instantiate "SpecialUnitaryGauge" object & its gauge field
  result = SpecialUnitaryGauge[n](actionPolicy: actionPolicy)
  result.newLatticeField(l,info)
  result.u = result.l[].newComplexGaugeLinks(n)

  # Set action coefficients
  case result.actionPolicy:
    of Wilson: result.actionParams = GaugeActionCoeffs(plaq: beta)
    of Symanzik: result.actionParams = gaugeActRect(beta, C1Symanzik)
    of Iwasaki: result.actionParams = gaugeActRect(beta, C1Iwasaki)
    of DoublyBlockedWilson: 
      result.actionParams = gaugeActRect(beta, C1DoublyBlockedWilson)
    of Rectangle: 
      if not info.hasKey("rectangle-coefficient"): echo GaugeWarning1
      let rectFac = info.hasKey("rectangle-coefficient")
        of true: info["rectangle-coefficient"].getFloat()
        of false: C1Symanzik
      result.actionParams = gaugeActRect(beta, rectFac)
    of Adjoint:
      if not info.hasKey("adjoint-ratio"): echo GaugeWarning2
      let adjFac = case info.hasKey("adjoint-ratio")
        of true: info["adjoint-ratio"].getFloat()
        of false: BetaAOverBetaF
      result.actionParams = GaugeActionCoeffs(plaq: beta, adjplaq: beta*adjFac)

proc newU1Gauge*(l: Layout; info: JsonNode): auto = 
  ## Instantiates U1 gauge field
  ## Inputs:
  ##   l: lattice layout (src/layout/)
  ##   info: JsonNode (species information about U1 gauge field)
  ## Output:
  ##   result: SpecialUnitaryGauge object
  result = l.newSpecialUnitaryGauge(info,1)

proc newSU2Gauge*(l: Layout; info: JsonNode): auto = 
  ## Instantiates SU2 gauge field
  ## Inputs:
  ##   l: lattice layout (src/layout/)
  ##   info: JsonNode (species information about SU2 gauge field)
  ## Output:
  ##   result: SpecialUnitaryGauge object
  result = l.newSpecialUnitaryGauge(info,2)

proc newSU3Gauge*(l: Layout; info: JsonNode): auto = 
  ## Instantiates SU3 gauge field
  ## Inputs:
  ##   l: lattice layout (src/layout/)
  ##   info: JsonNode (species information about SU3 gauge field)
  ## Output:
  ##   result: SpecialUnitaryGauge object
  result = l.newSpecialUnitaryGauge(info,3)

proc newSU4Gauge*(l: Layout; info: JsonNode): auto = 
  ## Instantiates SU4 gauge field
  ## Inputs:
  ##   l: lattice layout (src/layout/)
  ##   info: JsonNode (species information about SU4 gauge field)
  ## Output:
  ##   result: SpecialUnitaryGauge object
  result = l.newSpecialUnitaryGauge(info,4)

proc newSUNGauge*(l: Layout; info: JsonNode; n: static[int]): auto = 
  ## Instantiates SUN gauge field
  ## Inputs:
  ##   l: lattice layout (src/layout/)
  ##   info: JsonNode (species information about SU4 gauge field)
  ##   n: number of colors
  ## Output:
  ##   result: SpecialUnitaryGauge object
  result = l.newSpecialUnitaryGauge(info,n)

proc action*(self: SpecialUnitaryGauge): float =
  ## Calculate gauge action for special unitary field
  ## Inputs:
  ##   self: SpecialUnitaryGauge object
  ## Outputs:
  ##   result (float): value of gauge action
  result = case self.actionPolicy
    of Adjoint: self.actionParams.actionA(self.u)
    else: self.actionParams.gaugeAction1(self.u)

proc force*[S](self: SpecialUnitaryGauge; f: seq[S]) =
  ## Calculate gauge force for special unitary field
  ## Inputs:
  ##   self: SpecialUnitaryGauge object
  ##   f: force field
  ## Warning: f is written over, not appended to
  case self.actionPolicy:
    of Adjoint: self.actionParams.forceA(self.u,f)
    else: self.actionParams.gaugeForce(self.u,f)