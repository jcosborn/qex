when defined(QioReader) or not defined(QioliteReader):
  static: echo "Using QIO reader"
  import readerQio
  export readerQio
else:
  static: echo "Using Qiolite reader"
  import readerQiolite
  export readerQiolite
