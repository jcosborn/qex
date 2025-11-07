import qex
import gauge/stoutsmear
import utilityprocs

type
  StoutLinks*[S] = ref object
    nstout*: int
    stout*: seq[StoutSmear[seq[S]]]
    su*: seq[seq[S]]
    sf*: seq[seq[S]]

proc newStoutLinks[S](
    l: Layout,
    nstout: int,
    rho: float,
    s: typedesc[S]
  ): StoutLinks[S] =
  result = StoutLinks[S](nstout: nstout)
  result.stout = newSeq[StoutSmear[seq[S]]]()
  result.su = newSeq[seq[S]]()
  result.sf = newSeq[seq[S]]()
  for idx in 0..<result.nstout+1: 
    if idx != result.nstout: result.stout.add l.newStoutSmear(rho)
    result.su.add l.newGauge()
    result.sf.add l.newGauge()

proc newStoutLinks*(
    l: Layout,
    nstout: int,
    rho: float,
  ): auto = l.newStoutLinks(nstout,rho,type(l.newGauge()[0]))

proc smear*(self: StoutLinks; u: auto) =
  self.su[0].setGauge(u)
  for idx in 0..<self.nstout:
    self.stout[idx].smear(self.su[idx],self.su[idx+1])

proc smearedAction1*(
    gc: GaugeActionCoeffs,
    stout: StoutLinks;
    u: auto
  ): float = 
  stout.smear(u)
  result = gc.gaugeAction1(stout.su[^1])

proc smearForce(
    self: StoutLinks; 
    gc: GaugeActionCoeffs; 
    f: auto
  ) =
  for idx in countdown(self.nstout,0):
    let xdi = self.nstout - idx
    if xdi == 0: gc.gaugeActionDeriv(self.su[^1],self.sf[^1])
    elif xdi == self.nstout: self.stout[idx].smearDeriv(f,self.sf[1])
    else: self.stout[idx].smearDeriv(self.sf[idx],self.sf[idx+1])

proc smearedGaugeForce*(
    gc: GaugeActionCoeffs;
    stout: StoutLinks;
    u,f: auto
  ) =
  stout.smear(u)
  stout.smearForce(gc,f)
  contractProjectTAH(u,f)


    