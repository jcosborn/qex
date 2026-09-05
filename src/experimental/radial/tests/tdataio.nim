#RUNCMD $RUN1

import std/[os, strutils, tables, tempfiles, unittest]
import ../meas/dataio

addOutputFormatter(newConsoleOutputFormatter(colorOutput = false))

suite "measurement TSV identity and publication":
  setup:
    let dir = createTempDir("radial_tdataio_", "")
    let meta = @[("format", "radial-meas-2"), ("mass", "0.125"), ("seed", "17"),
                 ("pure", "false"), ("seaNf", "2"), ("source", "checkpoint"),
                 ("traj", "4"), ("nnoise", "16"), ("npair", "8"),
                 ("flowTimes", "0,0.20000000000000001,0.40000000000000002"),
                 ("srcTimes", "@[0, 30]"), ("ward", "false"), ("disc", "false")]
  teardown:
    removeDir dir

  test "identity allows descriptive metadata and relocation":
    let path = dir/"cond.t4.tsv"
    writeTsv(path, meta & @[("ensemble", "/previous/location"), ("note", "extra")],
             ["x"], [@[1.0]])
    requireTsvMeta(path, meta)
    check outputsDone([path], meta)
    check readTsvMeta(path)["ensemble"] == "/previous/location"

  test "each changed producing field identifies its mismatch":
    let path = dir/"cond.t4.tsv"
    writeTsv(path, meta, ["x"], [@[1.0]])
    for i in 0..<meta.len:
      var want = meta
      want[i][1] &= "changed"
      var msg = ""
      try:
        requireTsvMeta(path, want)
      except ValueError as e:
        msg = e.msg
      check msg.contains(meta[i][0])
      check msg.contains(path)
      check msg.contains("-skipDone:false")

  test "old metadata is rejected with the missing field":
    let path = dir/"cond.t4.tsv"
    writeTsv(path, [("massConvention", "standard"), ("overlapRho", "1")],
             ["x"], [@[1.0]])
    var msg = ""
    try:
      requireTsvMeta(path, meta)
    except ValueError as e:
      msg = e.msg
    check msg.contains("format")
    check msg.contains("<missing>")

  test "single files and gluon and current groups require every member":
    for i, names in [@["cond"], @["scalars"], @["wspec"], @["currents"],
                      @["flow", "loops", "f2"], @["currents", "currdisc", "currtrace"]]:
      var paths: seq[string]
      for name in names: paths.add dir/($i & "_" & name & ".t4.tsv")
      check not outputsDone(paths, meta)
      for k, path in paths:
        writeTsv(path, meta, ["x"], [@[1.0]])
        check outputsDone(paths, meta) == (k == paths.high)
      for path in paths:
        removeFile path
        check not outputsDone(paths, meta)
        writeTsv(path, meta, ["x"], [@[1.0]])

  test "a missing member does not hide an incompatible companion":
    let paths = [dir/"flow.t4.tsv", dir/"loops.t4.tsv", dir/"f2.t4.tsv"]
    writeTsv(paths[1], [("format", "old")], ["x"], [@[1.0]])
    expect ValueError:
      discard outputsDone(paths, meta)

  test "invalid replacement leaves the published data intact":
    let path = dir/"cond.t4.tsv"
    writeTsv(path, meta, ["x"], [@[1.0]])
    let old = readFile(path)
    expect ValueError:
      writeTsv(path, meta, ["x", "y"], [@[2.0], @[3.0, 4.0]])
    check readFile(path) == old
    expect ValueError:
      writeTsv(path, meta, ["x"], [@[2.0], @[3.0]])
    check readFile(path) == old
    writeTsv(path, meta, ["x"], [@[2.0, 1.0/3.0]])
    check readTsv(path).cols == @[@[2.0, 1.0/3.0]]
    check outputsDone([path], meta)
    for tmp in walkFiles(dir/"*.tmp"): check false

  test "failed publication removes the temporary file":
    let path = dir/"occupied"
    createDir path
    writeFile(path/"keep", "previous data")
    expect OSError:
      writeTsv(path, meta, ["x"], [@[1.0]])
    check readFile(path/"keep") == "previous data"
    for tmp in walkFiles(dir/"*.tmp"): check false

  test "preformatted fields preserve unavailable values and text labels":
    let path = dir/"summary.tsv"
    writeTsv(path, meta, ["name", "value", "err"],
             [@["Delta_A", "Delta_F"], @["1.25", "-2"], @["nan", "0.125"]])
    let text = readFile(path)
    check text.endsWith("Delta_A\t1.25\tnan\nDelta_F\t-2\t0.125\n")
    check readTsvMeta(path)["format"] == "radial-meas-2"

  test "invalid text replacement preserves the published data":
    let path = dir/"summary.tsv"
    writeTsv(path, meta, ["value"], [@["nan"]])
    let old = readFile(path)
    for ch in ['\t', '\n', '\r']:
      expect ValueError:
        writeTsv(path, meta, ["value"], @[@["1" & ch & "2"]])
      check readFile(path) == old
    for tmp in walkFiles(dir/"*.tmp"): check false
