import os
import stdUtils

when existsEnv("QIODIR"):
  const qioDir = getEnv("QIODIR")
else:
  const homeDir = getHomeDir()
  const qioDir = homeDir & "lqcd/install/qio"
{. passC: "-I" & qioDir & "/include" .}
{. passL: "-L" & qioDir & "/lib -lqio -llime" .}
{. pragma:qio, header:"qio.h".}

const 
  DML_UNKNOWN* = - 1
  DML_SINGLEFILE* = 0
  DML_MULTIFILE* = 1
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

type
  QIO_String* {.qio.} = object
    string*:cstring
    length*:csize
    
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
  QIO_Reader* {.qio.} = object
  QIO_RecordInfo* {.qio.} = object
  
proc QIO_string_create*():ptr QIO_String {.qio.}
proc QIO_string_destroy*(qs:ptr QIO_String) {.qio.}

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
