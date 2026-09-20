## NNFT checkpoint and latent-angle I/O using synthetic temporary arrays.
import base/globals
setDefaultNc(1)
setVLENmax(4)

import qex
import io/arrays
import ../io as ftio
import std/[json, math, os, tempfiles, unittest]

proc values[T: SomeFloat](leaf,n: int): seq[T] =
  result = newSeq[T](n)
  for j in 0..<n:
    result[j] = T(float32((37*j+19*leaf) mod 1024-512)/1024'f32)

proc spec(file,dtype: string; shape: openArray[int]): JsonNode =
  %*{"file":file,"dtype":dtype,"shape": @shape,"order":"C","byte_order":"little"}

proc saveManifest(dir: string; man: JsonNode) =
  writeFile(dir/"manifest.json",pretty(man))

proc makeCheckpoint(dir: string): JsonNode =
  const shapes = [@[12],@[12,6,3,3],@[12],@[12,12,3,3],@[12,1,1]]
  result = %*{"version":1,"arrays":{}}
  for i in 0..<40:
    let arr = spec("param_" & $i & ".bin","float32",shapes[i mod 5])
    arrays.writeArray(dir,arr,values[float32](i,arraySize(arr)))
    result["arrays"]["param_" & $i] = arr
  saveManifest(dir,result)

proc run() =
  suite "NNFT application I/O":
    test "checkpoint order, shapes and precision promotion":
      let dir = createTempDir("qex-nnft-io-","")
      defer: removeDir(dir)
      discard makeCheckpoint(dir)
      let single = loadNnft(dir)
      let double = loadNnft(dir,float64)
      check single.len == 8 and double.len == 8
      for s in 0..<8:
        check single[s].layers.len == 2
        check single[s].layers[0].cin == 6 and single[s].layers[0].cout == 12
        check single[s].layers[1].cin == 12 and single[s].layers[1].cout == 12
        check single[s].layers[0].bias == values[float32](5*s,12)
        check single[s].layers[0].weights == values[float32](5*s+1,12*6*9)
        check single[s].layers[1].bias == values[float32](5*s+2,12)
        check single[s].layers[1].weights == values[float32](5*s+3,12*12*9)
        check single[s].scale == values[float32](5*s+4,12)
        check double[s].layers[0].bias == values[float64](5*s,12)
        check double[s].layers[0].weights == values[float64](5*s+1,12*6*9)
        check double[s].layers[1].bias == values[float64](5*s+2,12)
        check double[s].layers[1].weights == values[float64](5*s+3,12*12*9)
        check double[s].scale == values[float64](5*s+4,12)
        for k in 0..<9:
          check single[s].layers[0].offsets[k] == @[int32(k div 3-1),int32(k mod 3-1)]
          check double[s].layers[1].offsets[k] == single[s].layers[0].offsets[k]

    test "leaf shapes define the stage count, channel widths and kernels":
      let dir = createTempDir("qex-nnft-io-","")
      defer: removeDir(dir)
      # Two stages of three convolutions: 6 -> 8 (5x5) -> 8 (3x3) -> 12 (3x1).
      const shapes = [@[8],@[8,6,5,5],@[8],@[8,8,3,3],@[12],@[12,8,3,1],@[12,1,1]]
      var man = %*{"version":1,"arrays":{}}
      for i in 0..<14:
        let arr = spec("param_" & $i & ".bin","float32",shapes[i mod 7])
        arrays.writeArray(dir,arr,values[float32](i,arraySize(arr)))
        man["arrays"]["param_" & $i] = arr
      saveManifest(dir,man)
      let p = loadNnft(dir)
      check p.len == 2
      for s in 0..1:
        check p[s].layers.len == 3
        check p[s].layers[0].kernel == @[5,5]
        check p[s].layers[2].kernel == @[3,1]
        check p[s].layers[1].cin == 8 and p[s].layers[2].cout == 12
        check p[s].layers[1].weights == values[float32](7*s+3,8*8*9)
        check p[s].scale == values[float32](7*s+6,12)
      man["arrays"].delete("param_13")
      saveManifest(dir,man)
      expect ValueError: discard loadNnft(dir)

    test "checkpoint shape, storage type, version and byte count are checked":
      let dir = createTempDir("qex-nnft-io-","")
      defer: removeDir(dir)
      let man = makeCheckpoint(dir)
      let arr = man["arrays"]["param_1"]
      arr["shape"] = %(@[6,12,3,3])
      saveManifest(dir,man)
      expect ValueError: discard loadNnft(dir)
      arr["shape"] = %(@[12,6,3,3])
      arr["dtype"] = %"float64"
      arrays.writeArray(dir,arr,values[float64](1,arraySize(arr)))
      saveManifest(dir,man)
      expect ValueError: discard loadNnft(dir,float64)
      arr["dtype"] = %"float32"
      arrays.writeArray(dir,arr,values[float32](1,arraySize(arr)))
      man["version"] = %2
      saveManifest(dir,man)
      expect ValueError: discard loadNnft(dir)
      man["version"] = %1
      saveManifest(dir,man)
      writeFile(dir/arr["file"].getStr,"short")
      expect ValueError: discard loadNnft(dir)

    test "latent angles map global direction, row and column in both precisions":
      let dir = createTempDir("qex-nnft-io-","")
      defer: removeDir(dir)
      let lo = newLayout(@[8,16])
      let g = lo.newGauge
      for dtype in ["float32","float64"]:
        let arr = spec("latent.bin",dtype,@[2,8,16])
        var data = newSeq[float64](arraySize(arr))
        for j in 0..<data.len:
          data[j] = 0.43*sin(0.31*float(j))+0.07*cos(0.19*float(j))
        arrays.writeArray(dir,arr,data)
        # The file holds the dtype-rounded values; fast-math may fold an
        # in-process float32 round trip away.
        let stored = arrays.readArray[float64](dir,arr)
        let man = %*{"version":1,"arrays":{"latent":arr}}
        saveManifest(dir,man)
        loadLatent(g,dir,"latent")
        for d in 0..1:
          for s in 0..<lo.nSites:
            let r = lo.coords[0][s].int
            let c = lo.coords[1][s].int
            let a = stored[(d*8+r)*16+c]
            var re,im: float64
            re := g[d]{s}[0,0].re
            im := g[d]{s}[0,0].im
            check abs(re-cos(a)) < 2e-16
            check abs(im-sin(a)) < 2e-16
        loadAngles(g,dir,arr)

    test "latent errors preserve destination storage":
      let dir = createTempDir("qex-nnft-io-","")
      defer: removeDir(dir)
      let lo = newLayout(@[8,16])
      let g = lo.newGauge
      threads:
        for f in g: f := 1.0
      let arr = spec("latent.bin","float64",@[2,8,16])
      arrays.writeArray(dir,arr,values[float64](0,arraySize(arr)))
      let man = %*{"version":1,"arrays":{"latent":arr}}
      saveManifest(dir,man)
      expect ValueError: loadLatent(g,dir,"missing")
      arr["shape"] = %(@[2,16,8])
      expect ValueError: loadAngles(g,dir,arr)
      arr["shape"] = %(@[2,8,16])
      writeFile(dir/"latent.bin","short")
      expect ValueError: loadAngles(g,dir,arr)
      for d in 0..1:
        for s in 0..<lo.nSites:
          var re,im: float64
          re := g[d]{s}[0,0].re
          im := g[d]{s}[0,0].im
          check re == 1.0
          check im == 0.0

qexInit()
run()
qexFinalize()
