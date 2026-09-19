## Array/field mapping tests use generated data and temporary files.
import qex
import io/[arrays, arrayfields]
import std/[json, os, unittest]

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))
qexInit()
let scratch = getTempDir()/"qex-arrayfields-test"/($getCurrentProcessId())
createDir(scratch)

proc spec(name, dtype: string, shape: seq[int]): JsonNode =
  %*{"file":name,"dtype":dtype,"shape":shape,"byte_order":"little","order":"C"}

# Lattices lay out up to VLEN 16: four dimensions with even outer extents.
const dims = @[4,4,8,8]

proc run[T: SomeFloat]() =
  let dtype = when T is float32: "float32" else: "float64"
  let dir = scratch/dtype
  let lo = newLayout(dims)
  proc field(): auto =
    when T is float32: lo.RealS()
    else: lo.RealD()
  let fs = @[field(),field()]

  suite "real field array IO " & dtype:
    test "channel arrays load at global coordinates":
      let meta = spec("load.bin",dtype,@[2] & dims)
      var data: seq[T]
      for c in 0..1:
        for x in 0..<dims[0]:
          for y in 0..<dims[1]:
            for z in 0..<dims[2]:
              for t in 0..<dims[3]: data.add T(100*c+10*x+y)+T(z)/T(8)-T(t)/T(4)
      writeArray(dir,meta,data)
      loadFields(fs,dir,meta)
      for c in 0..1:
        for s in lo.sites:
          var value: T
          value := fs[c]{s}
          check value == T(100*c+10*lo.coords[0][s]+lo.coords[1][s])+T(lo.coords[2][s])/T(8)-T(lo.coords[3][s])/T(4)

    test "saving preserves channel order and global C-order coordinates":
      # One directory for every rank: the master creates and broadcasts it. The
      # readback relies on the barrier inside saveFields, not on one here.
      var shared = ""
      if lo.comm.isMaster:
        shared = getTempDir()/"qex-arrayfields-shared"/($getCurrentProcessId())
        createDir(shared)
      lo.comm.broadcast(shared)
      let meta = spec("save.bin",dtype,@[2] & dims)
      for c in 0..1:
        for s in lo.sites:
          fs[c]{s} := T(100*c+lo.coords[0][s])+T(lo.coords[1][s])/T(8)-T(lo.coords[2][s])/T(16)+T(lo.coords[3][s])/T(32)
      saveFields(fs,shared,meta)
      let back = @[field(),field()]
      loadFields(back,shared,meta)
      for c in 0..1:
        for s in lo.sites:
          var a, b: T
          a := fs[c]{s}
          b := back[c]{s}
          check a == b
      if lo.comm.isMaster:
        let data = readArray[T](shared,meta)
        var n = 0
        for c in 0..1:
          for x in 0..<dims[0]:
            for y in 0..<dims[1]:
              for z in 0..<dims[2]:
                for t in 0..<dims[3]:
                  check data[n] == T(100*c+x)+T(y)/T(8)-T(z)/T(16)+T(t)/T(32)
                  inc n
        check n == data.len
      lo.comm.barrier
      if lo.comm.isMaster: removeDir(shared)

    test "shape and field layout mismatches fail before reading or writing":
      expect ValueError: loadFields(fs,dir,spec("absent.bin",dtype,@[2,4,8,4,8]))
      expect ValueError: saveFields(fs,dir,spec("absent.bin",dtype,@[1] & dims))
      # Every rank validates the descriptor before the collective transfer.
      expect ValueError: saveFields(fs,dir,spec("absent.bin","int32",@[2] & dims))
      let other = newLayout(dims)
      let f = when T is float32: other.RealS() else: other.RealD()
      expect ValueError: saveFields(@[fs[0],f],dir,spec("absent.bin",dtype,@[2] & dims))

suite "mask array IO":
  test "masks map nonzero bytes to selected scalar sites":
    let lo = newLayout(@[4,4,8,8])
    let meta = spec("mask.bin","uint8",@[4,4,8,8])
    var data: seq[uint8]
    for x in 0..<4:
      for y in 0..<4:
        for z in 0..<8:
          for t in 0..<8: data.add (if x mod 3 == 1 and y mod 2 == 0 and z != t: 255'u8 else: 0'u8)
    writeArray(scratch,meta,data)
    let mask = readMask(lo,scratch,meta)
    for s in lo.sites:
      var value: float32
      value := mask{s}
      let want = lo.coords[0][s] mod 3 == 1 and lo.coords[1][s] mod 2 == 0 and lo.coords[2][s] != lo.coords[3][s]
      check value == (if want: 1'f32 else: 0'f32)
    expect ValueError: discard readMask(lo,scratch,spec("absent.bin","uint8",@[4,8,4,8]))

  test "mask IO selects one site of the lattice":
    let lo = newLayout(dims)
    let meta = spec("mask3.bin","uint8",dims)
    var data: seq[uint8]
    for x in 0..<dims[0]:
      for y in 0..<dims[1]:
        for z in 0..<dims[2]:
          for t in 0..<dims[3]: data.add (if x == 2 and y == 3 and z == 7 and t == 1: 2'u8 else: 0'u8)
    writeArray(scratch,meta,data)
    let mask = readMask(lo,scratch,meta)
    for s in lo.sites:
      var value: float32
      value := mask{s}
      let want = lo.coords[0][s] == 2 and lo.coords[1][s] == 3 and lo.coords[2][s] == 7 and lo.coords[3][s] == 1
      check value == (if want: 1'f32 else: 0'f32)

run[float32]()
run[float64]()
removeDir(scratch)
qexFinalize()
