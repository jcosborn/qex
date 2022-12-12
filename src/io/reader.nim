import os

when defined(QioReader) or not defined(QioliteReader):
  import readerQio
  export readerQio
else:
  import readerQiolite
  export readerQiolite
