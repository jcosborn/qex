when defined(QioWriter) or not defined(QioliteWriter):
  static: echo "Using QIO writer"
  import writerQio
  export writerQio
else:
  static: echo "Using Qiolite writer"
  import writerQiolite
  export writerQiolite
