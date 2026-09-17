import io/arrays
import std/[json, os, tempfiles, unittest]

addOutputFormatter(newConsoleOutputFormatter(colorOutput=false))

proc spec(dtype: string; shape: openArray[int]; file = "array.bin"): JsonNode =
  %*{"file":file,"dtype":dtype,"shape": @shape,"byte_order":"little","order":"C"}

proc run() =
  let dir = createTempDir("qex-arrays-","")
  defer: removeDir(dir)

  suite "JSON binary arrays":
    test "known float32 bytes decode exactly in both precisions":
      let s = spec("float32",[3])
      writeFile(dir/"array.bin","\x00\x00\x80\x3f\x00\x00\x20\xc0\x00\x00\x00\x3e")
      check readArray[float32](dir,s) == @[1'f32,-2.5'f32,0.125'f32]
      check readArray[float64](dir,s) == @[1.0,-2.5,0.125]

    test "known float64 bytes decode exactly and convert to requested precision":
      let s = spec("float64",[2])
      writeFile(dir/"array.bin","\x00\x00\x40\x00\x00\x00\xf0\x3f\x00\x00\x00\x00\x00\x00\x04\xc0")
      check readArray[float64](dir,s) == @[1.0+1.0/1073741824.0,-2.5]
      check readArray[float32](dir,s) == @[1'f32,-2.5'f32]

    test "writes use the specified little-endian dtype":
      let a = spec("float32",[2],"nested/single.bin")
      writeArray(dir,a,[1.0,-2.5])
      check readFile(dir/"nested/single.bin") == "\x00\x00\x80\x3f\x00\x00\x20\xc0"
      let b = spec("float64",[2],"nested/double.bin")
      writeArray(dir,b,[1'f32,-2.5'f32])
      check readFile(dir/"nested/double.bin") == "\x00\x00\x00\x00\x00\x00\xf0\x3f\x00\x00\x00\x00\x00\x00\x04\xc0"
      writeArray(dir,spec("float32",[2],"nested/skipped.bin"),[1.0,-2.5],write=false)
      check not fileExists(dir/"nested/skipped.bin")
      expect ValueError: writeArray(dir,a,[1.0],write=false)

    test "float arrays round trip through either file precision":
      let a = @[0.37'f32,-12.125'f32,1.0'f32/7.0'f32,0'f32]
      let b = @[0.370000000001,-12.125,1.0/7.0,0.0]
      for dtype in ["float32","float64"]:
        let s = spec(dtype,[2,2])
        writeArray(dir,s,a)
        check readArray[float32](dir,s) == a
        writeArray(dir,s,b)
        let actual = readArray[float64](dir,s)
        # Compare with the stored values: fast-math may fold an in-process
        # float32 round trip away.
        let stored = readArray[float32](dir,s)
        for i in 0..<b.len:
          check actual[i] == (if dtype == "float32": float64(stored[i]) else: b[i])

    test "uint8 masks preserve all eight bits":
      let s = spec("uint8",[2,2])
      writeFile(dir/"array.bin","\x00\x01\x80\xff")
      check readArray[uint8](dir,s) == @[0'u8,1'u8,128'u8,255'u8]
      writeArray(dir,s,[0'u8,3'u8,127'u8,255'u8])
      check readFile(dir/"array.bin") == "\x00\x03\x7f\xff"
      check readArray[uint8](dir,s) == @[0'u8,3'u8,127'u8,255'u8]

    test "scalar shapes and empty arrays have explicit element counts":
      let scalar = spec("float64",[])
      check arrayShape(scalar) == newSeq[int]()
      check arraySize(scalar) == 1
      writeArray(dir,scalar,[0.125])
      check readArray[float64](dir,scalar) == @[0.125]
      for dtype in ["float32","float64"]:
        let empty = spec(dtype,[3,0,7])
        check arrayShape(empty) == @[3,0,7]
        check arraySize(empty) == 0
        writeArray(dir,empty,newSeq[float64]())
        check readFile(dir/"array.bin").len == 0
        check readArray[float32](dir,empty).len == 0
        check readArray[float64](dir,empty).len == 0
      let empty = spec("uint8",[0])
      writeArray(dir,empty,newSeq[uint8]())
      check readArray[uint8](dir,empty).len == 0

    test "malformed shape and overflowing counts are rejected":
      expect ValueError: discard arrayShape(%*{"shape":[1,-1]})
      expect ValueError: discard arrayShape(%*{"shape":[1,2.0]})
      expect ValueError: discard arrayShape(%*{"shape":"2,3"})
      expect ValueError: discard arrayShape(%*{})
      expect ValueError: discard arrayShape(newJNull())
      expect ValueError: discard arraySize(spec("float32",[high(int),2]))
      expect ValueError: discard readArray[float32](dir,spec("float32",[high(int)]))
      check arraySize(spec("float32",[high(int),2,0])) == 0

    test "malformed encoding metadata is rejected":
      for key in ["file","dtype","byte_order","order"]:
        let s = spec("float32",[1])
        s.delete(key)
        expect ValueError: discard readArray[float32](dir,s)
        s[key] = %13
        expect ValueError: discard readArray[float32](dir,s)
      for change in [("file",""),("dtype","int32"),("byte_order","big"),("order","F")]:
        let s = spec("float32",[1])
        s[change[0]] = %change[1]
        expect ValueError: discard readArray[float32](dir,s)
        expect ValueError: writeArray(dir,s,[1'f32])

    test "array file names stay inside the directory":
      let inner = dir/"inner"
      createDir(inner)
      writeFile(dir/"sentinel.bin","sentinel")
      for name in ["../sentinel.bin","nested/../sentinel.bin",dir/"sentinel.bin"]:
        let s = spec("float32",[2],name)
        expect ValueError: discard readArray[float32](inner,s)
        expect ValueError: writeArray(inner,s,[1'f32,2'f32])
        expect ValueError: writeArray(inner,s,[1'f32,2'f32],write=false)
      check readFile(dir/"sentinel.bin") == "sentinel"
      let dots = spec("float32",[2],"nested/a..b.bin")
      writeArray(inner,dots,[1'f32,2'f32])
      check readArray[float32](inner,dots) == @[1'f32,2'f32]

    test "mask and floating dtypes cannot be interchanged":
      let mask = spec("uint8",[1])
      let real = spec("float32",[1])
      expect ValueError: discard readArray[float32](dir,mask)
      expect ValueError: discard readArray[uint8](dir,real)
      expect ValueError: writeArray(dir,mask,[1'f32])
      expect ValueError: writeArray(dir,real,[1'u8])

    test "short and extra payload bytes are rejected":
      let s = spec("float32",[2])
      writeFile(dir/"array.bin",newString(7))
      expect ValueError: discard readArray[float32](dir,s)
      writeFile(dir/"array.bin",newString(9))
      expect ValueError: discard readArray[float32](dir,s)
      writeFile(dir/"array.bin","sentinel")
      expect ValueError: writeArray(dir,s,[1'f32])
      check readFile(dir/"array.bin") == "sentinel"

    test "missing array files report an IO error":
      expect IOError: discard readArray[float32](dir,spec("float32",[1],"missing.bin"))

    test "manifest reading validates only the generic object boundary":
      writeFile(dir/"manifest.json","{\"schema\":\"application-owned\",\"arrays\":{}}")
      check readManifest(dir)["schema"].getStr == "application-owned"
      writeFile(dir/"manifest.json","[]")
      expect ValueError: discard readManifest(dir)
      writeFile(dir/"manifest.json","{")
      expect ValueError: discard readManifest(dir)
      removeFile(dir/"manifest.json")
      expect IOError: discard readManifest(dir)

run()
