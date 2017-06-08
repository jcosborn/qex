when defined(lapackLib):
  import eigens/hisqev
  export hisqev
  when defined(primmeDir):
    import eigens/qexPrimme
    export qexPrimme
