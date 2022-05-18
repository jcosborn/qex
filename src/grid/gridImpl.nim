import grid/GridDefs
export GridDefs

import grid/Init
export Init

import grid/GridUtils
export GridUtils

import grid/gridStag
export gridStag

import base/qexInternal

proc gridInit() =
  var argc {.importc: "cmdCount", global.}: cint
  var argv {.importc: "cmdLine", global.}: cstringArray
  Grid_init(argc.addr, argv.addr)
  qexSetFinalizeComms(false)

proc gridFini() =
  Grid_finalize()

qexGlobalPreInit.add gridInit
qexGlobalPostFinal.add gridFini
