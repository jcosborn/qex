#[
Author: Curtis Taylor Peterson

Contact: curtistaylorpetersonwork@gmail.com

Source file: mcmc/fields/gauge/gaugeDefinitions.nim

Description: 
  Defines basic objects/types for gauge fields.

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
import ../latticeFields
import ../../mcmc/mcmcErrorHandling

export qex
export gauge
export latticeFields
export mcmcErrorHandling

const
  BetaAOverBetaF* = -1.0/4.0
  C1Symanzik* = -1.0/12.0
  C1Iwasaki* = -0.331
  C1DoublyBlockedWilson* = -1.4088

type
  GaugeActionPolicy* = enum 
    Wilson, 
    Adjoint, 
    Rectangle, 
    Symanzik, 
    Iwasaki, 
    DoublyBlockedWilson

const
  GaugeError1* = """
  |------------------------- QEX error -------------------------|
    Gauge action not specified for gauge field(s). Specify in 
    JSON input for gauge field as "action": "<option>" with...
  
    <option> = 
    wilson,
    Wilson,
    rectangle,
    Rectangle,
    adjoint,
    Adjoint,
    symanzik,
    Symanzik,
    luescher-weiss,
    Luescher-Weiss,
    luscher-weiss,
    Luscher-Weiss,
    iwasaki,
    Iwasaki,
    doubly-blocked-wilson,
    doubly-blocked-Wilson,
    Doubly-Blocked-Wilson,
    dbw,
    dbw2

                          - rectangle -
    If you specify "rectangle" or "Rectangle", you must provide 
    the rectangle coefficient in JSON input as 
    "rectangle-coefficient": <rectangle-coefficient>. Otherwise,
    rectangle with default to "symanzik", "Symanzik", 
    "luescher-weiss", "Luescher-Weiss", "luscher-weiss", 
    "Luscher-Weiss" value of -1/12. 
    
                        - adjoint-plaquette -
    If you specify "adjoint" or "Adjoint", you must provide 
    beta_F/beta_A (beta_F = beta) as 
    "adjoint-ratio": <adjoint-ratio>. Otherwise, beta_F/beta_A 
    will default to -1/4. 
  |------------------------- QEX error -------------------------|
  """
  GaugeError2* = """
  |------------------------- QEX error -------------------------|
    Bare gauge coupling (beta) not specified for gauge 
    field(s). Specify in JSON input as "beta": <beta>.
  |------------------------- QEX error -------------------------|
  """
  GaugeError3* = """
  |------------------------- QEX error -------------------------|
    Invalid choice for gauge action. Possible choices are:

    wilson,
    Wilson,
    rectangle,
    Rectangle,
    adjoint,
    Adjoint,
    symanzik,
    Symanzik,
    luescher-weiss,
    Luescher-Weiss,
    luscher-weiss,
    Luscher-Weiss,
    iwasaki,
    Iwasaki,
    doubly-blocked-wilson,
    doubly-blocked-Wilson,
    Doubly-Blocked-Wilson,
    dbw,
    dbw2

                          - rectangle -
    If you specify "rectangle" or "Rectangle", you must provide 
    the rectangle coefficient in JSON input as 
    "rectangle-coefficient": <rectangle-coefficient>. Otherwise,
    rectangle with default to "symanzik", "Symanzik", 
    "luescher-weiss", "Luescher-Weiss", "luscher-weiss", 
    "Luscher-Weiss" value of -1/12. 
    
                        - adjoint-plaquette -
    If you specify "adjoint" or "Adjoint", you must provide 
    beta_F/beta_A (beta_F = beta) as 
    "adjoint-ratio": <adjoint-ratio>. Otherwise, beta_F/beta_A 
    will default to -1/4. 
  |------------------------- QEX error -------------------------|
  """

const
  GaugeWarning1* = """
  |------------------------ QEX warning ------------------------|
    Rectangle coefficient not specified for gauge field(s). 
    Defaulting to "Symanzik" (Luescher-Weiss) value of -1/12. 
    If you wish to simulate with a different rectangle factor, 
    specify "rectangle-coefficient": <rectangle-coefficient>
    in JSON input.
  |------------------------ QEX warning ------------------------|
  """
  GaugeWarning2* = """
  |------------------------ QEX warning ------------------------|
    Adjoint ratio beta_A/beta_F not specified for gauge 
    field(s). Defaulting to "Adjoint" value of -1/4. If you 
    wish to simulate with a different adjoint ratio, specify 
    "adjoint-ratio": <adjoint-ratio> in JSON input.
  |------------------------ QEX warning ------------------------|
  """

if isMainModule:
  echo GaugeError1
  echo GaugeError2
  echo GaugeError3
  echo GaugeWarning1
  echo GaugeWarning2