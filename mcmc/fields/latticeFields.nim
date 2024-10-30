#[
Author: Curtis Taylor Peterson

Contact: curtistaylorpetersonwork@gmail.com

Source file: mcmc/fields/latticeFields.nim

Description: 
  Defines primitive data types for fields of various gauge groups/representations, 
  along with their constructors. All such data types are used by other object
  for whatever manipulations are needed. If your application needs a certain 
  field data structure, this is the place to define it. 

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
import layout
import json

type
  DComplexMatrixV*[N:static[int]] = Field[VLEN,Color[MatrixArray[N,N,DComplexV]]]
  DComplexOneIdxRep*[N:static[int]] = Color[VectorArray[N,DComplexV]]
  DComplexBosonOneIdxV*[N:static[int]] = Field[VLEN,DComplexOneIdxRep[N]]
  DComplexSpinOneIdx*[M,N:static[int]] = Spin[VectorArray[M,DComplexOneIdxRep[N]]]
  DComplexFermionOneIdxV*[M,N:static[int]] = Field[VLEN,DComplexSpinOneIdx[M,N]]

type
  LatticeField*[L:static[int]] {.inheritable.} = object
    ## Lattice field object
    ## Attributes:
    ##   l: reference to lattice layout
    ##   info: JSON data type containing lattice field info
    l*: ref Layout[L]
    info*: JsonNode

proc newComplexGaugeLinks*(l: Layout; n: static[int]): auto =
  ## Creates gauge link field Umu(x) ϵ Mat(nxn,C)
  ## Inputs:
  ##   l: lattice layout (src/layout/)
  ##   n: number of "colors" for fundamental,
  ##      dimension of representation otherwise
  ## Ouputs:
  ##   result: gauge field (DComplexMatrixV)
  result = newSeq[DComplexMatrixV[n]](l.nDim)
  for mu in 0..<l.nDim: result[mu].new(l)

proc newComplexBosonOneIdxRep*(l: Layout; n: static[int]): auto = 
  ## Creates complex representation of gauge group phi(x) ϵ Vec(n,C)
  ## Inputs:
  ##   l: lattice layout (src/layout/)
  ##   n: dimension of one-index representation
  ## Ouputs:
  ##   result: boson field (DComplexBosonOneIdxV)
  result = DComplexBosonOneIdxV[n].new(l)

proc newComplexFermionOneIdxRep*(l: Layout; m,n: static[int]): auto =
  ## Creates complex representation of gauge group psi(x) ϵ Vec(m,C) ⊗ Vec(n,C)
  ## Inputs:
  ##   l: lattice layout (src/layout/)
  ##   n: dimension of one-index representation
  ## Ouputs:
  ##   result: fermion field (DComplexFermionOneIdxV)
  result = DComplexFermionOneIdxV[m,n].new(l)

proc `:=`*[N:static[int]](v,u: seq[DComplexMatrixV[N]]) =
  ## Saves "u" field into "mu" field, both of type DComplexMatrixV
  ## Inputs:
  ##   u: Input gauge field
  ##   v: Gauge field to hold "u" values
  threads:
    for mu in 0..<u.len: v[mu] := u[mu] 

proc norm2*[N:static[int]](u: seq[DComplexMatrixV[N]]): float =
  ## Calculate norm of link-like field
  ## Inputs:
  ##   u: Input link field
  ## Outputs:
  ##   result: sum_mu(|u[mu]|^2)
  var u2: float
  threads:
    var tu2 = 0.0
    for mu in 0..<u.len: tu2 += u[mu].norm2
    threadBarrier()
    threadMaster: u2 = tu2
  result = u2

if isMainModule:
  qexInit()

  # Create lattice
  var
    lat = @[8,8,8,8]
    lo = newLayout(lat)

  # Create test fields
  var
    u1 = lo.newComplexGaugeLinks(1)
    su3 = lo.newComplexGaugeLinks(3)
    su3fbos = lo.newComplexBosonOneIdxRep(3)
    su3fferm = lo.newComplexFermionOneIdxRep(4,3)

  qexFinalize()