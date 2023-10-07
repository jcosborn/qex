import os

when defined(QioWriter) or not defined(QioliteWriter):
  static: echo "using QIO writer"
  import writerQio
  export writerQio
else:
  static: echo "using Qiolite writer"
  import writerQiolite
  export writerQiolite
