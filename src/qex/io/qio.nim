import os
import qex/base/stdUtils

when existsEnv("QMPDIR"):
  const qmpDir = getEnv("QMPDIR")
else:
  const qmpDir = getHomeDir() & "lqcd/install/qmp"
when existsEnv("QIODIR"):
  const qioDir = getEnv("QIODIR")
else:
  const qioDir = getHomeDir() & "lqcd/install/qio"
{. passC: "-I" & qioDir & "/include" .}
{. passL: "-L" & qioDir & "/lib -lqio -llime -L" & qmpDir & "/lib -lqmp" .}
{. pragma: qio, header: "qio.h" .}

const
  DML_UNKNOWN* = -1
  DML_SINGLEFILE* = 0.cint
  DML_MULTIFILE* = 1.cint
  DML_PARTFILE* = 2
  DML_PARTFILE_DIR* = 3
  DML_FIELD* = 0
  DML_GLOBAL* = 1
  DML_HYPER* = 2
  DML_SERIAL* = 0
  DML_PARALLEL* = 1

const
  QIO_UNKNOWN* = DML_UNKNOWN
  QIO_SINGLEFILE* = DML_SINGLEFILE
  QIO_MULTIFILE* = DML_MULTIFILE
  QIO_PARTFILE* = DML_PARTFILE
  QIO_PARTFILE_DIR* = DML_PARTFILE_DIR
  QIO_FIELD* = DML_FIELD
  QIO_GLOBAL* = DML_GLOBAL
  QIO_HYPER* = DML_HYPER
  QIO_SERIAL* = DML_SERIAL
  QIO_PARALLEL* = DML_PARALLEL
  QIO_ILDGNO* = 0
  QIO_ILDGLAT* = 1

var QIO_TRUNC*{.importc,qio.}: cint

type
  QIO_String* {.qio.} = object
    string*: cstring
    length*: csize

type
  QIO_Layout* {.qio.} = object
    node_number*: proc (coords: ptr cint): cint {.nimcall.}
    node_index*: proc (coords: ptr cint): cint {.nimcall.}
    get_coords*: proc (coords: ptr cint; node: cint; index: cint) {.nimcall.}
    num_sites*: proc (node: cint): cint {.nimcall.}
    latsize*: ptr cint
    latdim*: cint
    volume*: csize
    sites_on_node*: csize
    this_node*: cint
    number_of_nodes*: cint
  DML_io_node_t* = proc (a2: cint): cint {.nimcall.}
  DML_master_io_node_t* = proc (): cint {.nimcall.}
  QIO_Filesystem* {.qio.} = object
    number_io_nodes*: cint
    `type`*: cint
    my_io_node*: DML_io_node_t
    master_io_node*: DML_master_io_node_t
    io_node*: ptr cint
    node_path*: cstringArray
  QIO_Iflag* {.qio.} = object
    serpar*: cint
    volfmt*: cint
  QIO_Oflag* {.qio.} = object
    serpar*: cint
    mode*: cint
    ildgstyle*: cint
    ildgLFN*: ptr QIO_String
  QIO_Reader* {.qio.} = object
  QIO_Writer* {.qio.} = object
  QIO_RecordInfo* {.qio.} = object

proc QIO_string_create*: ptr QIO_String {.qio.}
proc QIO_string_destroy*(qs: ptr QIO_String) {.qio.}
proc QIO_string_set*(qs: ptr QIO_String; string: cstring) {.qio.}

proc QIO_open_read*(xml_file: ptr QIO_String; filename: cstring;
                    layout: ptr QIO_Layout; fs: ptr QIO_Filesystem;
                    iflag: ptr QIO_Iflag): ptr QIO_Reader {.qio.}
proc QIO_close_read*(r: ptr QIO_Reader): cint {.qio.}
proc QIO_next_record*(r: ptr QIO_Reader): cint {.qio.}
proc QIO_read_record_info*(qr:ptr QIO_Reader;
                           record_info:ptr QIO_RecordInfo;
                           xml_record:ptr QIO_String): cint {.qio.}
proc QIO_get_record_date*(record_info:ptr QIO_RecordInfo):cstring {.qio.}
proc QIO_get_datatype*(record_info:ptr QIO_RecordInfo):cstring {.qio.}
proc QIO_get_precision*(record_info:ptr QIO_RecordInfo):cstring {.qio.}
proc QIO_get_colors*(record_info:ptr QIO_RecordInfo):cint {.qio.}
proc QIO_get_spins*(record_info:ptr QIO_RecordInfo):cint {.qio.}
proc QIO_get_typesize*(record_info:ptr QIO_RecordInfo):cint {.qio.}
proc QIO_get_datacount*(record_info:ptr QIO_RecordInfo):cint {.qio.}
proc QIO_read*(r:ptr QIO_Reader; record_info: ptr QIO_RecordInfo;
               xml_record: ptr QIO_String;
               put:proc(buf:cstring; index:csize; count:cint; arg:pointer)
               {.nimcall.};
               datum_size: csize; word_size: cint; arg: pointer):cint {.qio.}
proc QIO_get_reader_latdim*(`in`:ptr QIO_Reader):cint {.qio.}
proc QIO_get_reader_latsize*(`in`:ptr QIO_Reader):ptr cArray[cint] {.qio.}

proc QIO_open_write*(xml_file: ptr QIO_String; filename: cstring; volfmt: cint;
                     layout: ptr QIO_Layout; fs: ptr QIO_Filesystem;
                     oflag: ptr QIO_Oflag): ptr QIO_Writer {.qio.}
proc QIO_close_write*(`out`: ptr QIO_Writer): cint {.qio.}
