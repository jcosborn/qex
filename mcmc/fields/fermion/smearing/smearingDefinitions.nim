#[
Author: Curtis Taylor Peterson

Contact: curtistaylorpetersonwork@gmail.com

Source file: mcmc/fields/fermion/linkSmearing/smearingDefinitions.nim

Description: 
  Defines various data types and strings for generic smearing
  operations under mcmc/fields/fermion/smearing

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
import gauge
import gauge/[hypsmear,stoutsmear]
import ../../latticeFields
import ../../gauge/unitary/unitaryGauge
import ../../../mcmc/mcmcErrorHandling

export qex
export gauge
export hypsmear
export stoutsmear
export latticeFields
export unitaryGauge
export mcmcErrorHandling

type
  SmearingPolicy* = enum
    HypercubicSmearing,
    StoutSmearing,
    NoSmearing

const
  SmearingError1* = """
  |------------------------- QEX error -------------------------|
    Smearing option not chosen. Specify in JSON input for
    gauge smearing as "smearing": "<option>" with...

    <option> = 
    stout,
    Stout,
    nhyp
    nHYP,
    hypercubic,
    Hypercubic

                            - stout - 
    If "stout" or "Stout" is specified, then you must provide 
    "rho" coefficient as "smearing-coefficient" = <rho>; 
    otherwise, rho will default 0.1. You can also optionally 
    specify "smearing-levels" = <# smears> for multiple stout 
    smears. If you alteratively provide 
    "smearing-coefficients" = [<rho_1>,<rho_2>,...] with
    [<rho_1>,<rho_2>,...] an array of size <# smears>, each
    stout smear will have its own value for rho.

                            - nhyp - 
    If "nhyp", "nHYP", "hypercubic",or "Hypercubic" is 
    specified, you must provide alpha_1, alpha_2, and alpha_3 
    parameters as 
    "smearing-coefficients": [<alpha_1>,<alpha_2>,<alpha_3>];
    otherwise, they will default to alpha_1 = 0.4, 
    alpha_2 = 0.5, and alpha_3 = 0.5. Note that multiple levels
    of stout smearing is not currently supported, but could be
    if it is requested.
  |------------------------- QEX error -------------------------|
  """
  SmearingError2* = """
  |------------------------- QEX error -------------------------|
    Chosen smearing option is either invalid or unavailable.
    Please choose from the following:
    
    stout,
    Stout,
    nhyp
    nHYP,
    hypercubic,
    Hypercubic

                            - stout - 
    If "stout" or "Stout" is specified, then you must provide 
    "rho" coefficient as "smearing-coefficient" = <rho>; 
    otherwise, rho will default 0.1. You can also optionally 
    specify "smearing-levels" = <# smears> for multiple stout 
    smears. If you alteratively provide 
    "smearing-coefficients" = [<rho_1>,<rho_2>,...] with
    [<rho_1>,<rho_2>,...] an array of size <# smears>, each
    stout smear will have its own value for rho.

                            - nhyp - 
    If "nhyp", "nHYP", "hypercubic",or "Hypercubic" is 
    specified, you must provide alpha_1, alpha_2, and alpha_3 
    parameters as 
    "smearing-coefficients": [<alpha_1>,<alpha_2>,<alpha_3>];
    otherwise, they will default to alpha_1 = 0.4, 
    alpha_2 = 0.5, and alpha_3 = 0.5. Note that multiple levels
    of stout smearing is not currently supported, but could be
    if it is requested.
  |------------------------- QEX error -------------------------|
  """

proc convSeqToStr(sq: seq[float]): string = 
  result = "["
  for idx,sqv in sq: 
    if idx == 0: result = result & $(sqv)
    else: result = result & ", " & $(sqv)
  result = result & "]"

template instantiationLogMessage*(
    policy: SmearingPolicy;
    coeffs: seq[float];
    nsmears: int
  ): untyped =
  var 
    str = case policy 
      of HypercubicSmearing,StoutSmearing: 
        "Instantiated " & $(policy) & " with parameters = " & coeffs.convSeqToStr
      of NoSmearing: "Instantiated " & $(policy)
  str

proc identity*[S](f,chain:seq[S]) =
  discard

if isMainModule:
  echo SmearingError1
  echo SmearingError2