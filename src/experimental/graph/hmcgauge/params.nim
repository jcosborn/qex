import qex
from os import fileExists
import config
import integrator

type
  HmcGaugeInputs* = object
    lat*: seq[int]
    beta*: float
    seed*: uint
    config*: RunConfig

template installHmcGaugeParams*() =
  installStandardParams()
  echoParams()
  processHelpParam()

proc readHmcGaugeInputs*(): HmcGaugeInputs =
  letParam:
    gaugefile = ""
    savefile = "config"
    savefreq = 0
    lat =
      if fileExists(gaugefile):
        getFileLattice gaugefile
      else:
        if gaugefile.len > 0:
          qexWarn "Nonexistent gauge file: ", gaugefile
        @[8,8,8,16]
    beta = 5.4
    dt = 0.025
    trajsThermo = 0
    trajsTrain = 50
    trajsTrainlrWarm = 10
    trajsInfer = 0
    lrmax = 1.0
    lrmin = 0.0001
    weightDecay = 0.0
    seed:uint = 1234567891
    gintalg = "2MN"
    # Legacy parameter name; this is an integrator coefficient vector, not
    # always a single lambda coefficient.
    lambda = newSeq[float]()
    gsteps = 4
    alwaysAccept:bool = 0

  result.lat = lat
  result.beta = beta
  result.seed = seed
  let integratorKind = parseIntegratorKind(gintalg)
  let integratorCoeffValues = lambda
  result.config = RunConfig(
    gaugefile: gaugefile,
    savefile: savefile,
    dt: dt,
    lrmax: lrmax,
    lrmin: lrmin,
    weightDecay: weightDecay,
    trajsThermo: trajsThermo,
    trajsTrain: trajsTrain,
    trajsTrainlrWarm: trajsTrainlrWarm,
    trajsInfer: trajsInfer,
    savefreq: savefreq,
    gsteps: gsteps,
    integratorCoeffs: parseIntegratorCoeffs(integratorKind, integratorCoeffValues),
    alwaysAccept: alwaysAccept)

template persistHmcGaugeParams*() =
  processSaveParams()
  writeParamFile()
