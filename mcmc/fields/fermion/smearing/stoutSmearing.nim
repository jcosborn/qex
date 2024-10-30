#[
Author: Curtis Taylor Peterson

Contact: curtistaylorpetersonwork@gmail.com

Source file: mcmc/fields/fermion/linkSmearing/stoutSmearing.nim

Description: 
  Procedures for stout smearing based on src/gauge/stoutsmear.nim. See
  the following for more details:
    - src/gauge/stoutsmear.nim
    - "Analytic smearing of SU⁡(3) link variables in lattice QCD"
      - Morningstar, C., Peardon, M., Phys. Rev. D 69, 054501
      - DOI: https://doi.org/10.1103/PhysRevD.69.054501

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

type Stout*[N:static[int]] = ref object
  ## Stout smearing object
  ## Attributes:
  ##   wilsonAction: Wilson action coefficients for smearing;
  ##     see comment below
  ##   rho: Stout smearing "rho" parameter (su = exp(rho*Q)*u)
  ##   su: Smeared gauge link
  ## Comments:
  ##   - Stout smear is done with respect to Wilson "flow" 
  ##     action; see Eqn. 4.43 of doi:10.1007/978-94-024-0999-4
  wilsonAction: GaugeActionCoeffs
  rho: float
  su,sf,expasf: seq[DComplexMatrixV[N]]

proc newStout*[N](u: seq[DComplexMatrixV[N]]; rho: float): auto =
  ## Stout object constructor
  ## Inputs:
  ##   u: Gauge field (../../latticeFields.nim)
  ##   rho: Stout smear "rho" parameter (su = exp(rho*Q)*u)
  result = Stout[N](rho: rho, wilsonAction: GaugeActionCoeffs(plaq: 1.0))
  result.su = u[0].l.newComplexGaugeLinks(N)
  result.sf = u[0].l.newComplexGaugeLinks(N)
  result.expasf = u[0].l.newComplexGaugeLinks(N)

proc smear(self: Stout; u,q: auto; derivative: bool = false) =
  ## Smear gauge link "u" with one level of stout smearing
  ## Inputs:
  ##   Stout: Stout object
  ##   u: Gauge field to be smeared (../../latticeFields.nim)
  ##   q: Lie algebra element entering smearing as su = exp(rho*q)*u
  ## Comments:
  ##   - Based on src/gauge/stoutsmear.nim
  let alpha = -self.rho*u[0][0].nrows.float # "rho" rescaled by Nc
  var # su,sf,expasf are references to Stout field attributes
    su = self.su 
    sf = self.sf
    expasf = self.expasf
  self.wilsonAction.gaugeActionDeriv(u,q)
  threads:
    for mu in 0..<q.len:
      for s in q[mu]:
        let staple = u[mu][s]*q[mu][s].adj
        var temp {.noinit.}: evalType(q[mu][s])
        temp.projectTAH staple
        if derivative: sf[mu][s] := temp
        temp := exp(alpha*temp)
        if derivative: expasf[mu][s] := temp
        su[mu][s] := temp*u[mu][s]

proc derivative(self: Stout; derivative,chain,q,tu: auto) =
  ## Calculates derivative of one stout-smeared field w.r.t. another field
  ## Inputs:
  ##   self: Stout object
  ##   derivative: Field holding derivative
  ##   chain: Field derivative being calculated with respect to
  ##   q,tu: Convenience fields used for calculation
  ## Comments:
  ##   - Based on src/gauge/stoutsmear.nim
  ##   - derivative & chain must be distinct fields
  ##   - Conventions:
  ##     c† = c† exp(a sf)
  ##     D[exp(a sf)]† = d/d(exp(a sf))[exp(a sf) u] c†
  ##     D[exp(a sf)]† = u c†
  ##     D[a sf]† = d/d(a sf)[exp(a sf)] D[exp(a sf)]†
  ##     D[sf]† = d/d(f)[a sf] D[a sf]†
  ##     D[u]† = d/d(u)[exp(a sf) u]
  let alpha = -self.rho*derivative[0][0].nrows.float # "rho" rescaled by Nc
  let # f,exp are references to Stout field attributes
    su = self.su
    sf = self.sf
    expasf = self.expasf
  threads:
    for mu in 0..<q.len:
      for s in q[mu]: 
        var temp {.noinit.}: type(q[mu][s])
        temp := chain[mu][s]*su[mu][s].adj
        derivative[mu][s] := alpha*expDeriv(alpha*sf[mu][s],temp)
  gaugeForceDeriv(self.su,derivative,derivative,q,tu)
  threads:
    for mu in 0..<q.len:
      for s in q[mu]: derivative[mu][s] += expasf[mu][s].adj*chain[mu][s]

proc smear*[N:static[int]](
    stout: seq[Stout[N]]; 
    u,su: auto; 
    derivative: bool = false
  ) =
  ## Smear gauge link "u" with multiple levels of stout smearing
  ## Inputs:
  ##   stout: sequence of Stout objects
  ##   u: Gauge field to be smeared (../../latticeFields.nim)
  ##   su: Smeared gauge field
  var q = u[0].l.newComplexGaugeLinks(N)
  su := u
  for smear in stout:
    smear.smear(su,q,derivative=derivative)
    su := smear.su
  qexLog "Smeared links with StoutSmearing"

proc smearGetDerivative*[N:static[int],S](
    stout: seq[Stout[N]]; 
    u,su: seq[S]
  ): auto =
  ## Smears gauge links and returns proc for smearing force
  ## Inputs:
  ##   stout: Sequence of Stout objects
  ##   u: Unsmeared field
  ##   su: Smeared field
  stout.smear(u,su,derivative=true)
  proc smearedForce(f,chain:seq[S]) =
    ## Calculates stout-smeared force
    ## Inputs:
    ##   f: Smeared force
    ##   chain: Field derivative being calculate with respect to
    ## Comments:
    ##   - See comments under "derivative" for details
    ##   - If using to calculate dS/dU for some stout-smeared
    ##     action "S", make sure "chain" is output of 
    ##     already-calculated force with respect to the 
    ##     outermost smeared link
    ##   - Output is not rephased and not projected to its
    ##     traceless/anti-Hermitian part; this must be done after
    var
      (q,tu,sfa,sfb) = (
        f[0].l.newComplexGaugeLinks(N),
        f[0].l.newComplexGaugeLinks(N),
        f[0].l.newComplexGaugeLinks(N),
        f[0].l.newComplexGaugeLinks(N)
      )
    for idx in countdown(stout.len,0):
      let xdi = stout.len - idx
      if xdi == 0: sfa := chain
      elif xdi == stout.len: stout[idx].derivative(f,sfa,q,tu)
      else: 
        stout[idx].derivative(sfb,sfa,q,tu)
        sfa := sfb
    qexLog "Smeared force with StoutSmearing"
  result = smearedForce