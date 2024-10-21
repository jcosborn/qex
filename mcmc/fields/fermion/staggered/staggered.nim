"""
Author: Curtis Taylor Peterson

Contact: curtistaylorpetersonwork@gmail.com

Source file: mcmc/fields/fermion/staggered/staggered.nim

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
"""

import qex
import physics/[stagD,stagSolve]
import ../operators/diracOperators
import ../../latticeFields

export latticeFields

const
  ActionCGTol = 1e-20
  ForceCGTol = 1e-12
  ActionMaxCGIter = 10000
  ForceMaxCGIter = 10000

type
  StaggeredDiracOperator*[S] = ref object of DiracOperator[S]
    rephased*: bool
    shifter*: seq[Shifter[S,W]]

proc zeroOdd*[T](phi: T) =
  threads: phi.odd := 0

proc zero*[T](phi: T) =
  threads: phi := 0

proc zero*[T](phis: seq[T]) =
  threads: 
    for phi in phis: phi := 0