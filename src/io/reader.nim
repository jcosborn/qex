when defined(QioReader):
  static: echo "using QIO reader"
  import readerQio
  export readerQio
else:
  static: echo "Using Qiolite reader"
  import readerQiolite
  export readerQiolite
