when defined(lapackLib):
  import eigens/hisqev
  export hisqev
  when defined(primmeDir):
    import eigens/qexPrimme, eigens/qexPrimmeSvd
    export qexPrimme, qexPrimmeSvd
