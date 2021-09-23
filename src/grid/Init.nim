import GridDefs
## ************************************************************************************
##
##     Grid physics library, www.github.com/paboyle/Grid
##
##     Source file: ./lib/Init.h
##
##     Copyright (C) 2015
##
## Author: Peter Boyle <paboyle@ph.ed.ac.uk>
## Author: paboyle <paboyle@ph.ed.ac.uk>
##
##     This program is free software; you can redistribute it and/or modify
##     it under the terms of the GNU General Public License as published by
##     the Free Software Foundation; either version 2 of the License, or
##     (at your option) any later version.
##
##     This program is distributed in the hope that it will be useful,
##     but WITHOUT ANY WARRANTY; without even the implied warranty of
##     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##     GNU General Public License for more details.
##
##     You should have received a copy of the GNU General Public License along
##     with this program; if not, write to the Free Software Foundation, Inc.,
##     51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
##
##     See the full license in the file "LICENSE" in the top level distribution directory
## ***********************************************************************************
##   END LEGAL

## ignored statement

proc Grid_init*(argc: ptr cint; argv: ptr cstringArray) {.
    importcpp: "Grid::Grid_init(@)", header: "Grid/Grid.h".}
proc Grid_finalize*() {.importcpp: "Grid::Grid_finalize(@)", header: "Grid/Grid.h".}
##  internal, controled with --handle

proc Grid_sa_signal_handler*(sig: cint; si: ptr siginfo_t; `ptr`: pointer) {.
    importcpp: "Grid::Grid_sa_signal_handler(@)", header: "Grid/Grid.h".}
proc Grid_debug_handler_init*() {.importcpp: "Grid::Grid_debug_handler_init(@)",
                                header: "Grid/Grid.h".}
proc Grid_quiesce_nodes*() {.importcpp: "Grid::Grid_quiesce_nodes(@)",
                           header: "Grid/Grid.h".}
proc Grid_unquiesce_nodes*() {.importcpp: "Grid::Grid_unquiesce_nodes(@)",
                             header: "Grid/Grid.h".}
proc GridDefaultSimd*(dims: cint; nsimd: cint): Coordinate {.
    importcpp: "Grid::GridDefaultSimd(@)", header: "Grid/Grid.h".}
proc GridDefaultLatt*(): var Coordinate {.importcpp: "Grid::GridDefaultLatt(@)",
                                      header: "Grid/Grid.h".}
proc GridDefaultMpi*(): var Coordinate {.importcpp: "Grid::GridDefaultMpi(@)",
                                     header: "Grid/Grid.h".}
proc GridThreads*(): var cint {.importcpp: "Grid::GridThreads(@)",
                            header: "Grid/Grid.h".}
proc GridSetThreads*(t: cint) {.importcpp: "Grid::GridSetThreads(@)",
                             header: "Grid/Grid.h".}
proc GridLogTimestamp*(a1: cint) {.importcpp: "Grid::GridLogTimestamp(@)",
                                header: "Grid/Grid.h".}
proc GridLogLayout*() {.importcpp: "Grid::GridLogLayout(@)", header: "Grid/Grid.h".}
##  Common parsing chores

proc GridCmdOptionPayload*(begin: cstringArray; `end`: cstringArray; option: string): string {.
    importcpp: "Grid::GridCmdOptionPayload(@)", header: "Grid/Grid.h".}
proc GridCmdOptionExists*(begin: cstringArray; `end`: cstringArray; option: string): bool {.
    importcpp: "Grid::GridCmdOptionExists(@)", header: "Grid/Grid.h".}
proc GridCmdVectorIntToString*[VectorInt](vec: VectorInt): string {.
    importcpp: "Grid::GridCmdVectorIntToString(@)", header: "Grid/Grid.h".}
proc GridCmdOptionCSL*(str: string; vec: var stdvector[string]) {.
    importcpp: "Grid::GridCmdOptionCSL(@)", header: "Grid/Grid.h".}
proc GridCmdOptionIntVector*[VectorInt](str: string; vec: var VectorInt) {.
    importcpp: "Grid::GridCmdOptionIntVector(@)", header: "Grid/Grid.h".}
proc GridCmdOptionInt*(str: var string; val: var cint) {.
    importcpp: "Grid::GridCmdOptionInt(@)", header: "Grid/Grid.h".}
proc GridParseLayout*(argv: cstringArray; argc: cint; latt: var stdvector[cint];
                     simd: var stdvector[cint]; mpi: var stdvector[cint]) {.
    importcpp: "Grid::GridParseLayout(@)", header: "Grid/Grid.h".}
proc printHash*() {.importcpp: "Grid::printHash(@)", header: "Grid/Grid.h".}
## ignored statement
