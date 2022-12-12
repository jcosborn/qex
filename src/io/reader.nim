import os

when defined(QioReader):
  import readerQio
  export readerQio
else:
  import readerQiolite
  export readerQiolite
