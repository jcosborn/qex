#[
Author: Curtis Taylor Peterson

Contact: curtistaylorpetersonwork@gmail.com

Source file: mcmc/fields/fermion/linkSmearing/nhypSmearing.nim

Description: 
  Wrapper procedures for nHYP smearing (src/gauge/hypsmear). See the
  following for more details:
    - src/gauge/hypsmear.nim
    - "Hypercubic smeared links for dynamical fermions"
      - Hasenfratz, A., Hoffmann, R., Schaefer, S., JHEP05(2007)029
      - DOI: 10.1088/1126-6708/2007/05/029

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
import gauge/hypsmear

proc smear*[S](nhyp: var tuple[parameters:HypCoefs,info:PerfInfo]; u,v: seq[S]) =
  ## Smearing u links & store result in v links
  ## Inputs:
  ##   nhyp: Tuple of information for nHYP smearing (src/gauge/hypsmear.nim)
  ##   u: Input gauge field (unaffected)
  ##   v: Smeared u links
  ## Comments:
  ##   - Wraps smearing proc in src/gauge/hypsmear.nim
  nhyp.parameters.smear(u,v,nhyp.info)
  qexLog "Smeared links with HypercubicSmearing"

proc smearGetDerivative*[S](
    nhyp: var tuple[parameters:HypCoefs,info:PerfInfo]; 
    u,v: seq[S]
  ): auto =
  ## Smear u links, store result in v links, return smeared force proc
  ## Inputs:
  ##  nhyp: Tuple of information for nHYP smearing (src/gauge/hypsmear.nim)
  ##   u: Input gauge field (unaffected)
  ##   v: Smeared u links
  ## Output:
  ##   result: Proc taking in force and returning its smeared counterpart
  ## Comments:
  ##   - Wraps "smearGetForce" in src/gauge/hypsmear.nim
  let smearForce = nhyp.parameters.smearGetForce(u,v,nhyp.info)
  proc smearedForce(f,chain:seq[S]) =
    ## Smears force with nHYP smearing
    ## Inputs:
    ##   f: Smeared force
    ##   chain: Field derivative being calculate with respect to
    ## Comments:
    ##   - fₓₚₜ ← chainₘₖₕ d/dUₓₚₜ^*[Vₘₖₕ(U)^*] + chainₘₕₖ^* d/dUₓₚₜ^*[Vₘₕₖ(U)]
    ##   - If using to calculate dS/dU for some nHYP-smeared
    ##     action "S", make sure "chain" is output of 
    ##     already-calculated force with respect to the 
    ##     outermost smeared link
    ##   - Output is not rephased and not projected to its
    ##     traceless/anti-Hermitian part; this must be done after
    f.smearForce(chain)
    qexLog "Smeared force with HypercubicSmearing"
  qexLog "Smeared links with HypercubicSmearing"
  result = smearedForce
