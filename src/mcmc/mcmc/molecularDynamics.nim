import ../mcmcTypes
import ../utilities/gaugeUtils
import ../fields/gaugeFields
import ../fields/staggeredFields
import ../actions/latticeActionUtils

import gauge/[hypsmear, stoutsmear]

template smearForce(self: LatticeAction; f: auto): untyped = 
  self.smear.sf(f,f)

template smearForce(self: LatticeSubAction; f: auto): untyped = 
  self.smear[].sf(f,f)

proc getStep[A](self: A): seq[float] =
  result = newSeq[float]() 
  for sAction in self.subActions:
    if sAction.running: result.add sAction.absCumSum/sAction.absScale

proc setStep[A](self: A; minabsdt, scale: float): float =
  result = 0.0
  for sAction in self.subActions:
    case (sAction.running): #and (not sAction.cycleCompleted):
      of true:
        let abscdt = sAction.absCumSum/sAction.absScale
        if abs(minabsdt - abscdt) <= Small64: # Mod for force-grad.?
          let nMD = sAction.updates.len div 2
          sAction.timeStep = 2*(sAction.iStep).mod(nMD)
          sAction.spaceStep = sAction.timeStep + 1

          sAction.vdtau = sAction.updates[sAction.spaceStep].dtau*scale
          result = abscdt*scale
          if not sAction.solo: sAction.includeInStep = true

          sAction.iStep += 1
          sAction.timeStep = 2*(sAction.iStep).mod(nMD)
          sAction.spaceStep = sAction.timeStep + 1
          
          sAction.cumSum += sAction.updates[sAction.timeStep].dtau
          sAction.absCumSum += abs(sAction.updates[sAction.timeStep].dtau)
          if sAction.iStep == sAction.steps*nMD: sAction.running = false
        else: sAction.vdtau = 0.0
      of false: sAction.vdtau = 0.0

proc getNumRunning[A](self: A): int = 
  result = 0
  for sAction in self.subActions:
    if sAction.running: result += 1

proc smearGaugeForce(self: var Smearing) =
  if not self.smeared:
    case self.smearing:
      of Hypercubic: 
        self.sf = smearGetForce(self.nhyp, self.u[], self.su, self.nhypInfo)
      of Stout: discard
      of NoSmearing: discard
    self.setMatterBoundaryConditions
    self.smeared = true

template smearGaugeForce(self: LatticeSubAction): untyped =
  self.smear[].smearGaugeForce

proc rephaseForce(f: auto) =
  f.stagPhase
  threads:
    for mu in 0..<f.len:
      for s in f[mu].odd: f[mu][s] *= -1

proc projectTrclssAntiHrm(f: auto; u: auto) =
  threads:
    for mu in 0..<f.len:
      for s in f[mu]:
        var temp {.noinit.}: typeof(f[0][0])
        temp := f[mu][s]*u[mu][s].adj
        projectTAH(f[mu][s], temp)

proc fermionHeatbath(self: LatticeSubAction) =
  case self.pField.field:
    of GaugeField: discard
    of StaggeredMatterField:
      self.rephase
      self.pRNG[].randomComplexGaussian(self.D[][StaggeredMatterField].stagPsi)
      self.pField.getStaggeredField(self.D[][StaggeredMatterField])
    of WilsonMatterField: discard
    of DummyField: discard
  if not self.solo:
    for sAction in self.subActions: sAction.fermionHeatbath

proc fermionHeatbath*(self: LatticeAction) =
  for sAction in self.subActions: sAction.fermionHeatbath

proc zero(f: auto) =
  threads:
    for mu in 0..<f.len:
      for s in f[mu]: f[mu][s] := 0.0

proc updateMomentum(p: auto; f: auto; dtau: float) =
  threads:
    for mu in 0..<f.len:
      for s in f[mu]: p[mu][s] -= dtau*f[mu][s]

proc updateMomentum(p: auto; f: auto) =
  threads:
    for mu in 0..<f.len:
      for s in f[mu]: p[mu][s] -= f[mu][s]

proc getGaugeForce[A](self: A; f: auto; p: auto) =
  var 
    gaugeUpdated = false
    dtau = 0.0
  proc search(f: auto; self: LatticeSubAction) =
    case self.pField.field:
      of GaugeField: 
        if abs(self.vdtau) > Small64:
          gaugeForce(self.pField, f)
          gaugeUpdated = true
          dtau = self.vdtau
      else: discard
  for sAction in self.subActions: f.search(sAction)
  if gaugeUpdated: p.updateMomentum(f,dtau)


proc getStaggeredMatterForce(f: auto; self: LatticeSubAction): bool = 
  result = false
  case self.pField.field:
    of GaugeField,WilsonMatterField,DummyField: discard
    of StaggeredMatterField:
      if abs(self.vdtau) > Small64:
        self.smearGaugeForce
        self.rephase
        stagForce(self.pField, self.D[][StaggeredMatterField], self.vdtau, f)
        result = true

proc getWilsonMatterForce(f: auto; self: LatticeSubAction): bool = 
  result = false

proc getMatterForce(
    actions: seq[LatticeAction]; 
    u: auto; 
    f: auto; 
    p: auto
  ) =
  for action in actions:
    var 
      update = false
      updateStag = false
      updateWils = false
    case action.action:
      of PureGauge: discard
      of GaugeMatter,PureMatter:
        for sAction in action.subActions:
          if abs(sAction.vdtau) > Small64: update = true
        if update: f.zero

        for sAction in action.subActions:
          if f.getStaggeredMatterForce(sAction): updateStag = true
        if updateStag: f.rephaseForce
        
        for sAction in action.subActions:
          if f.getWilsonMatterForce(sAction): updateWils = true
        
        if (updateStag) or (updateWils):
          f.setMatterBoundaryConditions(action)
          action.smearForce(f)
          f.projectTrclssAntiHrm(u)
          p.updateMomentum(f)

proc getMatterForce(
    actions: seq[LatticeSubAction]; 
    u: auto; 
    f: auto; 
    p: auto
  ) =
  var 
    update = false
    updateStag = false
    updateWils = false

  for action in actions:
    case action.action:
      of PureGauge: discard
      of GaugeMatter,PureMatter:
        for sAction in action.subActions:
          if abs(sAction.vdtau) > Small64: update = true
  if update: f.zero

  for action in actions:
    case action.action:
      of PureGauge: discard
      of GaugeMatter,PureMatter:
        for sAction in action.subActions:
          if f.getStaggeredMatterForce(sAction): updateStag = true
  if updateStag: f.rephaseForce

  for action in actions:
    case action.action:
      of PureGauge: discard
      of GaugeMatter,PureMatter:
        for sAction in action.subActions:
          if f.getWilsonMatterForce(sAction): updateWils = true
  
  if (updateStag) or (updateWils): # Have a common smearing
    f.setMatterBoundaryConditions(actions[0])
    actions[0].smearForce(f)
    f.projectTrclssAntiHrm(u)
    p.updateMomentum(f)

proc mdStep[A](
    actions: seq[A]; 
    dtau: float; 
    u: auto; 
    f: auto; 
    p: auto
  ) =
  var 
    nested = false
    hasMatter = false

  # Gauge field update
  for action in actions:
    for sAction in action.subActions:
      if not sAction.solo: nested = true
  case nested:
    of true: 
      var subActions = newSeq[type(actions[0].subActions[0])]()
      for action in actions:
        for sAction in action.subActions:
          if not sAction.solo:
            if sAction.includeInStep: 
              subActions.add sAction
              sAction.includeInStep = false
      case subActions.len > 0:
        of true: subActions.trajectory(u,f,p,type(subActions[0]),dtau)
        of false: u.updateGauge(p,dtau)
    of false: u.updateGauge(p,dtau)

  # Momentum update
  for action in actions:
    case action.action:
      of PureGauge: discard
      of GaugeMatter, PureMatter:
        action.smear.rephased = false
        action.smear.smeared = false
        hasMatter = true
  for action in actions:
    case action.action:
      of PureGauge, GaugeMatter: action.getGaugeForce(f,p)
      of PureMatter: discard
  if hasMatter: actions.getMatterForce(u,f,p)

  #[
  var ts = @[dtau]
  for action in actions:
    for sAction in action.subActions:
      ts.add sAction.vdtau
  echo ts
  ]#


template trajectory*[A](
    actions: seq[A]; 
    u: auto;
    f: auto;
    p: auto;
    a: typedesc[A];
    scale: float
  ): untyped =
  var 
    step = 0
    nRunning = 0

    minabsdt = Large64

    pTime = 0.0
    cTime = 0.0

    pdtau = 0.0
    dtau = 0.0
    
    updateSpace = false

  while true:
    minabsdt = Large64
    nRunning = 0
    updateSpace = false

    for action in actions:
      let absdts = action.getStep
      for absdt in absdts:
        if absdt <= minabsdt: minabsdt = absdt
    for action in actions: cTime = action.setStep(minabsdt, scale)
    dtau = cTime - pTime + pdtau

    for action in actions: nRunning += action.getNumRunning

    for action in actions:
      for sAction in action.subActions:
        if abs(sAction.vdtau) > Small64: updateSpace = true

    if (updateSpace) or (nRunning == 0): actions.mdStep(dtau,u,f,p)

    if nRunning == 0:
      for action in actions:
        for sAction in action.subActions:
          sAction.running = true
          sAction.iStep = 0
          sAction.timeStep = 0
          sAction.spaceStep = 0
          sAction.cumSum = sAction.updates[0].dtau
          sAction.absCumSum = abs(sAction.cumSum)
      break
    else:
      pdtau = case updateSpace
        of true: 0.0
        of false: dtau
      pTime = cTime