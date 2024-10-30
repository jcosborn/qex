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

import ../gaugeDefinitions
import json

type
  SpecialUnitaryGauge*[L,N:static[int]] = ref object of LatticeField[L]
    ## SpecialUnitaryGauge object
    ## Attributes:
    ##   u: gauge field (mcmc/fields/latticeFields.nim)
    ##   actionPolicy: named gauge action (mcmc/fields/gauge/gaugeDefinitions.nim)
    ##   actionParams: gauge action coefficients (src/gauge/gaugeAction.nim)
    u*: seq[DComplexMatrixV[N]] 
    actionPolicy: GaugeActionPolicy 
    actionParams: GaugeActionCoeffs 

# Main "SpecialUnitaryGauge" constructor
proc newSpecialUnitaryGauge[L:static[int]](
    l: Layout[L]; 
    info: JsonNode;
    n: static[int];
  ): auto =
  ## SpecialUnitaryGauge object constructor
  ## Inputs:
  ##   l: lattice layout (src/layout/)
  ##   info: JsonNode (species information about U1 or SU(N) gauge field)
  ##   n: number of "colors"
  ## Output:
  ##   result: SpecialUnitaryGauge object
  
  var 
    actionPolicy: GaugeActionPolicy
    beta: float

  if not info.hasKey("action"): qexError GaugeError1
  case info["action"].getStr():
    of "wilson","Wilson": actionPolicy = Wilson
    of "rectangle","Rectangle": actionPolicy = Rectangle
    of "adjoint","Adjoint": actionPolicy = Adjoint
    of "symanzik","Symanzik": actionPolicy = Symanzik
    of "luescher-weiss","Luescher-Weiss": actionPolicy = Symanzik
    of "luscher-weiss", "Luscher-Weiss": actionPolicy = Symanzik
    of "iwasaki","Iwasaki": actionPolicy = Iwasaki
    of "doubly-blocked-wilson","doubly-blocked-Wilson","Doubly-Blocked-Wilson":
      actionPolicy = DoublyBlockedWilson
    of "dbw","dbw2": actionPolicy =  DoublyBlockedWilson
    else: 
      echo GaugeError3
      throwError info["action"].getStr() & " not a valid action"

  if not info.hasKey("beta"): qexError GaugeError2
  beta = info["beta"].getFloat()

  # Instantiate "SpecialUnitaryGauge" object & its gauge field
  result = SpecialUnitaryGauge[L,n](actionPolicy: actionPolicy, info: info)
  new(result.l); result.l[] = l;
  result.u = result.l[].newComplexGaugeLinks(n)

  # Set action coefficients
  case result.actionPolicy:
    of Wilson: result.actionParams = GaugeActionCoeffs(plaq: beta)
    of Rectangle: 
      if not info.hasKey("rectangle-coefficient"): echo GaugeWarning1
      let rectFac = case info.hasKey("rectangle-coefficient")
        of true: info["rectangle-coefficient"].getFloat()
        of false: C1Symanzik
      result.actionParams = gaugeActRect(beta, rectFac)
    of Adjoint:
      if not info.hasKey("adjoint-ratio"): echo GaugeWarning2
      let adjFac = case info.hasKey("adjoint-ratio")
        of true: info["adjoint-ratio"].getFloat()
        of false: BetaAOverBetaF
      result.actionParams = GaugeActionCoeffs(plaq: beta, adjplaq: beta*adjFac)
    else: # Nim compiler workaround
      if result.actionPolicy == Symanzik: 
        result.actionParams = gaugeActRect(beta, C1Symanzik)
      elif result.actionPolicy == Iwasaki:
        result.actionParams = gaugeActRect(beta, C1Iwasaki)
      elif result.actionPolicy == DoublyBlockedWilson:
        result.actionParams = gaugeActRect(beta, C1DoublyBlockedWilson)

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
  
if isMainModule:
  qexInit()

  # Create lattice
  var
    lat = @[8,8,8,8]
    lo = newLayout(lat)

  # Instantiate gauge field objects
  var
    info = %* {
      "action": "Wilson",
      "beta": 7.0
    }
    su2 = lo.newSU2Gauge(info)
    su3 = lo.newSU3Gauge(info)
    su2f = lo.newComplexGaugeLinks(2)
    su3f = lo.newComplexGaugeLinks(3)

  # Set u to unit gauge
  su2.u.unit
  su3.u.unit

  # Echo unit gauge action
  echo su2.action, " ", su3.action

  # Calculate unit gauge force
  su2.force(su2f)
  su3.force(su3f)

  qexFinalize()