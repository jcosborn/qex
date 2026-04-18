import nimibook

var book = initBookWithToc:
  entry("Introduction", "index", numbered = false)
  entry("Installation", "install")
  entry("Building", "build")
  section("QEX Tutorial", "qex/index.nim"):
    #entry("Getting Started", "getting_started.nim")
    entry("Simple example", "simple_example.nim")
    entry("Basic Operations", "basic_ops.nim")

nimibookCli(book)
