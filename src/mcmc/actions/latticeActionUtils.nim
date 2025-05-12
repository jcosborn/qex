import ../mcmcTypes
import ../fields/gaugeFields
import ../fields/staggeredFields

import gauge/[hypsmear, stoutsmear]

template su(self: LatticeAction): untyped = self.smear.su

template u(self: LatticeAction): untyped = self.smear.u[]

template nhyp(self: LatticeAction): untyped = self.smear.nhyp

template nhypInfo(self: LatticeAction): untyped = self.smear.nhypInfo

template su(self: LatticeSubAction): untyped = self.smear[].su

template u(self: LatticeSubAction): untyped = self.smear[].u[]

template nhyp(self: LatticeSubAction): untyped = self.smear[].nhyp

template nhypInfo(self: LatticeSubAction): untyped = self.smear[].nhypInfo

proc setMatterBoundaryConditions*(su: auto; bc: string) =
  threads:
    for mu in 0..<su.len:
      if $(bc[mu]) == $"a":
        tfor i, 0..<su[mu].l.nSites:
          if su[mu].l.coords[mu][i] == su[mu].l.physGeom[mu]-1:
            su[mu]{i} *= -1.0

proc setMatterBoundaryConditions*(self: Smearing) =
  self.su.setMatterBoundaryConditions(self.bc)

proc setMatterBoundaryConditions*(f: auto; self: LatticeAction) =
  let bc = self.smear.bc
  f.setMatterBoundaryConditions(bc)

proc setMatterBoundaryConditions*(f: auto; self: LatticeSubAction) =
  let bc = self.smear[].bc
  f.setMatterBoundaryConditions(bc)

template setMatterBoundaryConditions(self: LatticeAction): untyped =
  self.smear.setMatterBoundaryConditions

template setMatterBoundaryConditions(self: LatticeSubAction): untyped = 
  self.smear[].setMatterBoundaryConditions

proc rephase(self: Smearing) =
  if not self.rephased: 
    stagPhase(self.su)
    self.rephased = true

template rephase*(self: LatticeSubAction): untyped =
  self.smear[].rephase

proc resetRephase(self: Smearing) =
  self.rephased = false

template resetRephase*(self: LatticeAction): untyped =
  case self.action:
    of PureGauge: discard
    of GaugeMatter,PureMatter: resetRephase(self.smear)

proc smearGauge*[A](self: A) =
  case self.action:
    of PureGauge: discard
    of GaugeMatter,PureMatter:
      case self.smearing:
        of Hypercubic: smear(self.nhyp, self.u, self.su, self.nhypInfo)
        of Stout: discard
        of NoSmearing: discard
      self.setMatterBoundaryConditions

proc getAction(self: LatticeSubAction): float =
  case self.pField.field:
    of GaugeField: self.currentAction = gaugeAction(self.pField)
    of StaggeredMatterField:
      self.rephase
      self.currentAction = staggeredAction(self.pField, self.D[][StaggeredMatterField])
    of WilsonMatterField: discard
    of DummyField: discard
  result = self.currentAction
  if not self.solo:
    for sAction in self.subActions: result += sAction.getAction

proc getAction*(self: LatticeAction): float =
  result = 0.0
  for sAction in self.subActions: result += sAction.getAction
  self.resetRephase