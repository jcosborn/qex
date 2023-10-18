when defined(QioWriter):
  static: echo "using QIO writer"
  import writerQio
  export writerQio
else:
  static: echo "Using Qiolite writer"
  import writerQiolite
  export writerQiolite
