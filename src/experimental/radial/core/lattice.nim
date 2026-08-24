## S^2 x R lattice and the free-limit couplings.
## Normative reference: doc/02-formulation.md sections 3, 5, 6; doc/04-interfaces.md section 4.

import std/math
import geom
export geom

type
  Lat* = ref object
    sph*: Sphere
    nt*: int
    at*: float
    nsite*: int              ## sph.nv * nt
    kap*: seq[float]         ## kappa_e = 2 A_e/(abar l_e), A_e the exact kite area (IV.2);
                             ## the flat identity (l*_1 + l*_2)/abar holds only to O(abar^2)
    kapT*: seq[float]        ## kappa'_y = A_y/(abar at)                          (IV.2)
    volw*: seq[float]        ## A_y at, the local volume element                  (IV.11)
    volbar*: float           ## mean of volw

proc newLat*(sph: Sphere, nt: int, at: float): Lat =
  result = Lat(sph: sph, nt: nt, at: at, nsite: sph.nv*nt,
               kap: newSeq[float](sph.ne), kapT: newSeq[float](sph.nv),
               volw: newSeq[float](sph.nv))
  for e in 0..<sph.ne:
    result.kap[e] = 2.0*sph.edges[e].area/(sph.abar*sph.edges[e].len)
  var v = 0.0
  for y in 0..<sph.nv:
    result.kapT[y] = sph.area[y]/(sph.abar*at)
    result.volw[y] = sph.area[y]*at
    v += result.volw[y]
  result.volbar = v/float(sph.nv)

func sIdx*(l: Lat, v, t: int): int = v + l.sph.nv*floorMod(t, l.nt)
func eIdx*(l: Lat, e, t: int): int = e + l.sph.ne*floorMod(t, l.nt)
func tIdx*(l: Lat, v, t: int): int = v + l.sph.nv*floorMod(t, l.nt)

func betaFace*(l: Lat, f: int, g2: float): float =
  ## beta_tri = at/(g2 A_tri)   (IV.26)
  l.at/(g2*l.sph.faces[f].area)

func betaEdge*(l: Lat, e: int, g2: float): float =
  ## beta_l = 2 A_e/(g2 l_e^2 at)   (IV.26)
  2.0*l.sph.edges[e].area/(g2*l.sph.edges[e].len*l.sph.edges[e].len*l.at)

func asOverAt*(l: Lat): float =
  ## Doubler condition (section 6) requires abar/at >= 4/3.
  l.sph.abar/l.at

func maxM*(l: Lat): float =
  ## 0.9 M_0 with M_0 = min(4/sqrt(3), sqrt(3) abar/at)   (IV.10)
  0.9*min(4.0/sqrt(3.0), sqrt(3.0)*l.sph.abar/l.at)
