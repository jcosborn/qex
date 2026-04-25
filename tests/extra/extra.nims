# Extra Tests
#
# extraTest(sourceFile, runArgs="")
#
# If a directory with the same base name as the source with an extra
# 't' in front appears in this directory, then the test script will
# run the 'run' script within that directory, otherwise it will just
# run the executable.

extraTest("gauge/wflow.nim")
extraTest("examples/staghmc_sh.nim")
extraTest("examples/staghmc_shf.nim")
extraTest("examples/generichmc.nim")
