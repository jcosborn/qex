/**
 * @file HighlyImprovedStaggeredFermionImpl.h
 * @brief Interface for implementation of highly improved staggered fermions (HISQ)
 * @author Curtis Taylor Peterson
 * @details
 * This header file is meant to act as an interface for both Grid and the 
 * MILC codebase to utilize the "highly improved staggered quark" action
 * within Grid. The interface itself is composed of a single core component class
 * 
 * * HighlyImprovedStaggeredFermionImpl: Class providing support for highly
 *   improved staggered fermions in Grid. Exposes methods for fat7/asqtad 
 *   smearing and its derivative with optional inclusion of Lepage term needed 
 *   for second level of HISQ smearing. 
 * 
 * References:
 *   * Quantum EXpressions (QEX): https://github.com/jcosborn/qex
 *   * QOPQDP [SciDAC]: https://github.com/usqcd-software/qopqdp
 *   * SIMULATEeQCD: https://github.com/LatticeQCD/SIMULATeQCD
 *   * Follana, E. et al.: https://doi.org/10.1103/PhysRevD.75.054502
 *   * MILC Collaboration (2010): https://doi.org/10.1103/PhysRevD.82.074501
 * 
 * Acknowledgements:
 *   Curtis Taylor Peterson would like to thank James Osborn for developing/testing
 *   the implementations of HISQ in QEX and QOPQDP, from which the "fast" option for 
 *   the fat7/asqtad smearing and derivativative are based and have been tested against.
 *   
 *   This material is based upon work supported by the U.S. Department of Energy, 
 *   Office of Science, Office of Advanced Scientific Computing Research, Scientific 
 *   Discovery through Advanced Computing (SciDAC) program.
*/

#pragma once 

#ifndef QCD_UTILS_HISQ_IMPL_H
#define QCD_UTILS_HISQ_IMPL_H 

#include <Grid/Grid.h>
#include "PeriodicTransporters.h"
#include "UnitaryProjection.h"

NAMESPACE_BEGIN(Grid);

//
// macros
//

// loop excluding no indices
#define HISQLOOP0(exec) for (int mu = 0; mu < Nd; ++mu) {exec;} \

// loop excluding one index
#define HISQLOOP1(exec) for (int nu = 0; nu < Nd; ++nu) \
  {if (nu == mu) {continue;} else exec;}                \

// loop excluding two indices
#define HISQLOOP2(exec) for (int i = 0; i < Nd; ++i)   \
  {if ((i == mu) or (i == nu)) {continue;} else exec;} \

// loop excluding three indices
#define HISQLOOP3(exec) for (int j = 0; j < Nd; ++j)               \
  {if ((j == mu) or (j == nu) or (j == i)) {continue;} else exec;} \

// encapsulates conditional execution of Lepage calculations for cleanliness
#define HISQLEPAGE(exec) if (ctx.lepage != 0.0) { exec; } \

// encapsulates conditional execution of Naik calculations for cleanliness
#define HISQNAIK(exec) if (ctx.naik != 0.0) { exec; } \

// better notation for ternary expression
#define when(cond, valTrue, valFalse) ((cond) ? (valTrue) : (valFalse)) \

//
// useful consts
//

// usual 4D HISQ coefficients: defined for convenience
const double // Lepage & Naik
  LEPAGE = -1.0/8.0,
  NAIK = -1.0/24.0;
const double // fat-7
  F7L1 = -LEPAGE,
  F7L3 = -0.5*F7L1,
  F7L5 = -0.25*F7L3,
  F7L7 = 0.0625*NAIK;
const double // asqtad
  ASQL1 = -8.0*LEPAGE,
  ASQL3 = F7L3,
  ASQL5 = F7L5,
  ASQL7 = F7L7; 
const bool BACKUPSVD = true; // use backup SVD in unitary projection when applicable
const double // unitary projection parameters
  REUNITCUTOFF = 1e-20,           // base-level cutoff on eigenvalues
  REUNITDERIVCUTOFF = 5e-5,  // base-level cutoff on eigenvalues for derivative
  BACKUPSVDTOLERANCE = 1e-8;      // tolerance for triggering backup SVD

//
// convenient data structures
//

struct StagImplParams {
  Coordinate dirichlet; // Blocksize of dirichlet BCs
  int  partialDirichlet;
  AcceleratorVector<Complex, Nd> boundary_phases;
  StagImplParams()
  { partialDirichlet = 0; dirichlet.resize(0); boundary_phases.resize(Nd, 1.0); };
  StagImplParams(std::vector<Complex> phi): boundary_phases(phi)
  { partialDirichlet = 0; dirichlet.resize(0); };
};

struct HISFContext {
  /**
   * @brief Context for highly improved staggered fermion smearing and projection
   * @author Curtis Taylor Peterson
   * @details
   * Context structure for passing parameters to HISF smearing and projection methods.
   */
  // fat7 smearing parameters
  RealD c0, c1, c2, c3;

  // asqtad parameters
  RealD lepage;
  RealD naik;

  // unitary projection parameters
  bool backupSVD;
  RealD svdTolerance;
  RealD eigenCutoff;

  HISFContext(
    RealD c0, 
    RealD c1, 
    RealD c2, 
    RealD c3,
    RealD lepage, 
    RealD naik,
    bool backupSVD,
    RealD svdTolerance,
    RealD eigenCutoff
  ):
    c0(c0), 
    c1(c1), 
    c2(c2), 
    c3(c3),
    lepage(lepage), 
    naik(naik),
    backupSVD(backupSVD),
    svdTolerance(svdTolerance),
    eigenCutoff(eigenCutoff) { };
  
  HISFContext(RealD c0, RealD c1, RealD c2, RealD c3, RealD lepage, RealD naik):
    c0(c0), c1(c1), c2(c2), c3(c3), lepage(lepage), naik(naik) { };
  
  HISFContext(RealD c0, RealD c1, RealD c2, RealD c3): 
    c0(c0), c1(c1), c2(c2), c3(c3), lepage(0.0), naik(0.0) { };
  
  HISFContext(bool backupSVD, RealD svdTolerance, RealD eigenCutoff):
    backupSVD(backupSVD), svdTolerance(svdTolerance), eigenCutoff(eigenCutoff) { };
};

//
// highly improved staggered fermion implementation
//

template <class Gimpl>
class HighlyImprovedStaggeredFermionImpl: Gimpl {
/**
 * @brief Highly improved staggered fermion implementation in Grid
 * @author Curtis Taylor Peterson
 * @details
 * Highly improved staggered fermion (HISF/HISQ) implementation in Grid. Exposes 
 * methods for FNAL/MILC and non-FNAL/MILC implementations of Highly Improved Staggered
 * Fermions. This class boasts an implementation of HISQ that is low in communication 
 * overhead whilst retaining a low memory footprint by Grid's PaddedCell and 
 * GeneralLocalStencil through the PeriodicTransporters class.
 */
  
public: INHERIT_GIMPL_TYPES(Gimpl);

private:
  typedef typename Simd::scalar_type GridScalar;
  typedef typename std::vector<GaugeLinkField> GaugeLorentzField;
  typedef typename std::vector<ComplexField> StaggeredPhases;

public:
  GridBase* grid;
  GridBase* longGrid;

  PaddedCell cell;
  PaddedCell longCell;

  StagImplParams params;
  StaggeredPhases stagPhases;

private:
  void init(
    GridCartesian* grid,
    PaddedCell& cell, 
    PaddedCell& longCell, 
    bool calculateStaggeredPhases
  ) {
    assert(Nc == 3 && "HISQ only suppored for SU(3)");
    assert(Nd == 4 && "HISQ only supported for 4 dimensions");

    this->grid = (GridBase*)cell.grids.back();
    this->longGrid = (GridBase*)longCell.grids.back();

    if (calculateStaggeredPhases) calcStagPhases(stagPhases);
  }

public:
  HighlyImprovedStaggeredFermionImpl(
    GridCartesian* grid, 
    const StagImplParams params,
    bool calculateStaggeredPhases = true
  ):
    cell(1, grid), 
    longCell(2, grid),
    stagPhases(Nd, grid), 
    params(params)
  { init(grid, cell, longCell, calculateStaggeredPhases); }

  HighlyImprovedStaggeredFermionImpl(
    GridCartesian* grid,
    bool calculateStaggeredPhases = true
  ):
    cell(1, grid), 
    longCell(2, grid),
    stagPhases(Nd, grid), 
    params(StagImplParams({1, 1, 1, -1}))
  { init(grid, cell, longCell, calculateStaggeredPhases); }


private:
  void calcStagPhases(StaggeredPhases& phases) {
    /**
     * @brief HISQ gauge configuration constructor
     * @author Curtis Taylor Peterson, Peter Boyle
     * @details
     * Staggered phi follow the "MILC convention", which treats the fourth 
     * direction as the "time" coordinate:
     * (2a) eta_0 = (-1)^{x3}       <---| 
     * (2b) eta_1 = (-1)^{x3+x0}        | convention in
     * (2c) eta_2 = (-1)^{x3+x0+x1}     | this code
     * (2d) eta_3 = 1               <---|
     * Though awkard, this convention follows that of most texts on
     * relativity. This is opposed to the convention that one often finds in 
     * lattice gauge theory textbooks, where the "time" direction is x0:
     * (3a) eta_0 = 1
     * (3b) eta_1 = (-1)^{x0}
     * (3c) eta_2 = (-1)^{x0+x1}
     * (3d) eta_3 = (-1)^{x0+x1+x2}
     * Dirichlet boundary conditions are imposed by "rephasing" the links 
     * on the boundary, with "1" for periodic and "-1" for anti-periodic.
     */
    GridBase *grid = phases[0].Grid();
    Lattice<iScalar<vInteger>> x(grid), y(grid), t(grid);
    Lattice<iScalar<vInteger>> tx(grid), txy(grid), xyzt(grid), coor(grid);
    ComplexField phi(grid);

    LatticeCoordinate(x, 0);
    LatticeCoordinate(y, 1);
    LatticeCoordinate(t, 3);
    tx = t + x;
    txy = tx + y;

    for (int mu = 0; mu < Nd; mu++) { 
      // for boundary phi
      int N = grid->GlobalDimensions()[mu] - 1;
      auto bpha = params.boundary_phases[mu];
      GridScalar dirichlet(real(bpha), imag(bpha));

      // staggered phases x boundary phases
      LatticeCoordinate(coor, mu);
      phi = 1.0;
      if (mu == 0){phi = where(mod(t, 2) == (Integer)0,   phi, -phi);}
      if (mu == 1){phi = where(mod(tx, 2) == (Integer)0,  phi, -phi);}
      if (mu == 2){phi = where(mod(txy, 2) == (Integer)0, phi, -phi);}
      phases[mu] = where(coor == (Integer)N, dirichlet*phi, phi);
    }
  }

private:
  // wrap PokeIndex (insert gauge link field into gauge field layout)
  inline void insertLink(GaugeField& Uout, const GaugeLinkField& U, int mu)
  { PokeIndex<LorentzIndex>(Uout, U, mu); }

  // wrap PeekIndex (extract gauge link field form gauge field layout)
  inline const GaugeLinkField toLink(const GaugeField& U, int mu)
  { return PeekIndex<LorentzIndex>(U, mu); }

  // break up memory layout of gauge field into std::vector of gauge link fields
  inline const GaugeLorentzField toLorentz(const GaugeField& U) 
  { GaugeLorentzField u(Nd, grid); HISQLOOP0(u[mu] = toLink(U, mu)) return u; }

  // aggregate std::vector of gauge link fields into layout of gauge field
  inline const GaugeField toGauge(GaugeLorentzField& u) 
  { GaugeField U(grid); HISQLOOP0(insertLink(U, u[mu], mu)) return U; }

public:
  /** @brief rephases gauge links with staggered and Dirichlet boundary phases */
  void rephase(GaugeField& X, const GaugeField& W) 
  { for (int mu = 0; mu < Nd; ++mu) insertLink(X, stagPhases[mu]*toLink(W, mu), mu); }

public:
  void smear(
    GaugeField& X,
    GaugeField& WWW,
    const GaugeField& W,
    const HISFContext ctx
  ) {
    /**
     * @brief "Effective 3-link" fat7/asqtad + Lepage (optional) smear
     * @author Curtis Taylor Peterson
     * @details
     * Fat7/Asqtad links can be constructed recursively from 3-link smears. This 
     * insight is attributed to James Osborn (Argonne National Laboratory),
     * who used it to implement the Fat7/Asqtad smearing in QOPQDP and Qauntum
     * EXpressions (QEX). This method follows a similar pattern to James' 
     * implementation and combines it with Grid's padded cell. Note that many other 
     * codes, such as MILC, use a "path-based" approach by default, which tend to
     * less efficiently reuse computations from lower levels in the smearing.
     * 
     * Define
     * (1) Sν(U;n) := Uν(n) U(n+ν) Uν(n+μ)† + Uν(n-ν)† U(n-ν) U(n-ν+μ).
     * for some generic SU(N) field U(n). I call the first contribution the "upper 
     * staple" and the second the "lower staple". From (1), recursively define
     * (2) Sν(n) := Sν(Uμ;n),
     * (3) Sνj(n) := Sj(Sν;n),
     * and
     * (4) Sνij(n) := Si(Sνj;n).
     * Then an fat7 link Vμ(n) in four dimensions can be expressed as a composition 
     * of sums
     * (5) Vμ(n) = c0 Uμ(n)            [5a]
     *     + ∑_{ν≠μ}(c1 Sν(n)          [5b]
     *     + ∑_{i≠μ,ν}[c2 Sνi(n)       [5c]
     *     + ∑_{j≠μ,ν,i} c3 Sνij(n)])  [5d]
     * And that's it. The key thing to note is that we can reuse the staples
     * at lower levels to construct staples at higher levels.
     * 
     * The optional Lepage term adds an additional contribution to (5b):
     * (6) + lepage * Sν(Sν(n)),
     * where Sν(n) is treated as the "link" field in the Lepage staple. We can reuse
     * the code that calculates the symmetric staple, so long as we correct for the
     * additional 1-link terms that doing so adds to the smearing by compensating
     * in the 1-link coefficient. 
     * 
     * References:
     * * Quantum EXpressions (QEX): https://github.com/jcosborn/qex
     * * QOPQDP [SciDAC]: https://github.com/usqcd-software/qopqdp
     * * Follana, E. et al.: https://doi.org/10.1103/PhysRevD.75.054502
     */
    //PeriodicTransporters<Gimpl> w(when(ctx.naik == 0.0, cell, longCell), W);
    PeriodicTransporters<Gimpl> w(cell, W);
    GaugeLorentzField x(Nd, grid);
    GaugeLinkField s3(grid), s5(grid);
    
    HISQLOOP0(                                               // Eqn 5a
      x[mu] = (ctx.c0 - 6.0*ctx.lepage)*w.link(mu);
      HISQLOOP1(                                             // Eqn 5b
        s3 = w.staple(w.link(mu), mu, nu);
        x[mu] += ctx.c1*s3; 
        HISQLEPAGE(x[mu] += ctx.lepage*w.staple(s3, mu, nu)) // Eqn 6
        HISQLOOP2(                                           // Eqn 5c
          s5 = w.staple(s3, mu, i);
          x[mu] += ctx.c2*s5;
          HISQLOOP3(x[mu] += ctx.c3*w.staple(s5, mu, j))     // Eqn 5d
    ) ) )
    X = w.toTightGrid(toGauge(x));

    HISQNAIK(
      HISQLOOP0(x[mu] = w[mu].CovShift(w[mu].CovShift(w.link(mu), FORWARD), FORWARD))
      WWW = ctx.naik*w.toTightGrid(toGauge(x));
    )
  }

  void smear(GaugeField& X, GaugeField& WWW, const GaugeField& W) {
    HISFContext asqtadCtx(ASQL1, ASQL3, ASQL5, ASQL7, LEPAGE, NAIK); 
    smear(X, WWW, W, asqtadCtx); 
  }

  void smear(GaugeField& X, const GaugeField& W, const HISFContext ctx) 
  { smear(X, X, W, ctx); }

  void smear(GaugeField& X, const GaugeField& W) { 
    HISFContext fat7Ctx(F7L1, F7L3, F7L5, F7L7);  
    smear(X, X, W, fat7Ctx); 
  }

  void project(GaugeField& V, const GaugeField& U, const HISFContext ctx) { 
    UnitaryProjection<Gimpl> projection(
      ctx.eigenCutoff, 
      ctx.backupSVD, 
      ctx.svdTolerance
    );
    projection.project(V, U);
  }

  void project(GaugeField& V, const GaugeField& U, bool forDerivative = false) {
    UnitaryProjection<Gimpl> projection(
      when(forDerivative, REUNITDERIVCUTOFF, REUNITCUTOFF), 
      BACKUPSVD, 
      BACKUPSVDTOLERANCE
    );
    projection.project(V, U);
  }

public:
  void smearDerivative(
    GaugeField& dXdU,
    const GaugeField& dXdW,
    const GaugeField& dXdWWW,
    const GaugeField& W,
    const HISFContext ctx
  ) {
    /**
     * @brief "Effective 3-link" fat7/asqtad + Lepage (optional) smear
     * @author Curtis Taylor Peterson
     * @details
     * Calculates Wirtinger derivative of fat7/asqtad + Lepage links, including chain
     * rule. See the comments under `smear` method for details on fat7/asqtad + 
     * Lepage smearing. This method was originally implemented by translating the
     * "effective 3-link" force in QOP/QDP and Quantum EXpressions (QEX) written by 
     * James Osborn (Argonne National Laboratory) into Grid. It was then reworked to 
     * make full use of thevpadded cells in Grid, which resulted in an algorithm that 
     * is similar in spirit to that in QOP/QDP and QEX, but which deviates in how it 
     * reasons through each contribution to the derivative. The result is an 
     * algorithm that boasts a significant reduction in memory usage compared to 
     * QOP/QDP and QEX. This code has been tested against its translated counterpart,
     * which itself was thoroughly tested against QEX.
     * 
     * The full force from some smeared action S(Uμ) is
     * (1) δS/δUμ(n) = LieAlgebraProjection[ Uμ(n) ∂S/∂Uμ(n) ]
     * This method calculates ∂S/∂Uμ(n) (i.e., not δS/δUμ(n)); it is to be 
     * interpreted as a matrix Wirtinger derivative, with the rules
     * (2) [∂S/∂Uμ(n)]ij = ∂S/∂Uμ(n)ji,                      [2a]
     *     ∂Uμ(n)ij/∂Uν(m)kl = δ(i,k) δ(j,l) δ(μ,ν) δ(m,n),  [2b]
     *     ∂Uμ(n)/∂Uν†(m) = 0.                               [2c]
     * Taking Vμ(n) to be some fat7/asqtad-smeared link [Eqn. 5 in `smear` method], 
     * the chain rule gives
     * (2) ∂S/∂Uμ(n) = ∑_{δ,m} ∂S/∂Vδ(m) ∂Vδ(m)/∂Uμ(n)
     *               + ∑_{δ,m} ∂S/∂Vδ†(m) ∂Vδ†(m)/∂Uμ(n).
     * For simplicity, consider the top contribution to a single symmetric Sν(n): 
     * ν      ---🠢
     * 🠡      🠡   🠣   (Fig 1)
     * -🠢 μ   ....
     * The forward derivative of Sν(n) with respect to Uμ(n) will receive 
     * contributions from the leftmost link and the top link. To understand what
     * we're doing here, we need to consider these contributions "from the perspective
     * of the force". In other words, before we take the derivative, we think of the 
     * the contribution from the chain rule ∂S/∂Vδ as completing a plaquette 
     * (just with the like ∂S/∂Vδ pointing in the "wrong" direction); in Fig. 1, this 
     * is the dotted line below the staple. We then think of the contribution from 
     * the derivative on the left and top link by "rotating" the perspective in 
     * Fig. 1, such that the link that the derivative ∂/∂Uμ(n) hits lies on the 
     * x-axis where the contribution from ∂S/∂Vδ was. We then get
     * (3)            Uμ(n)
     *        🠠----   ---🠢 
     * ν      ⮾   🠡 + 🠡   🠣 (Fig 2)
     * 🠡      ----🠢    -⮾-
     * -🠢 μ   Uμ(n)
     *         [3a]    [3b]
     * where an "⮾" indicates a replacement with the chain rule; only when the 
     * center link is hit is the chain rule and the derivative aligned along 
     * the same direction. Additionally, notice that the orientation of the staples
     * that comprise the staple derivative are opposite to those in the origanl 
     * link; this is due to rotating the perspective to that of the force. To make 
     * diagrammatics work out without writing redundant procedures for the opposite
     * orientation, we work with the adjoint of the chain rule, calculate the 
     * derivative as if we were completing properly-oriented plaquettes, then correct
     * by returning the adjoint of the full derivative. The full contribution from the 
     * 3-link staples follows by repeating this procedure for the adjoint, reflection, 
     * and adjoint of the reflection. The 5- and 7-link staples can be constructed 
     * from the 3-link staples in a similar manner.
     * 
     * As the derivative is quite complicated, some of this may seem a bit vague.
     * Please feel free to contact Curtis Peterson (see email above) for any 
     * questions that may arise.
     * 
     * References:
     * * Quantum EXpressions (QEX): https://github.com/jcosborn/qex
     * * QOPQDP [SciDAC]: https://github.com/usqcd-software/qopqdp
     * * Follana, E. et al.: https://doi.org/10.1103/PhysRevD.75.054502
     */
    PeriodicTransporters<Gimpl> w(cell, W);
    GaugeLorentzField dxdw = toLorentz(w.toPaddedGrid(dXdW));
    GaugeLorentzField dxdu(Nd, grid);
    GaugeLinkField cnu(grid), ci(grid);
    GaugeLinkField snu(grid), si(grid), sj(grid);
    GaugeLinkField dsnu(grid), dsi(grid);

    HISQLOOP0( // fat7 + lepage (lepage won't execute if lepage == 0.0)
      dxdu[mu] = (ctx.c0 - 6.0*ctx.lepage)*dxdw[mu];
      HISQLOOP1(
        snu = ctx.c1*w.link(nu);
        dsnu = Zero();
        cnu = ctx.c1*dxdw[mu];
        HISQLEPAGE(
          snu += ctx.lepage*w.staple(w.link(nu), nu, mu);
          dsnu += ctx.lepage*w.staple(dxdw[nu], nu, mu);
          cnu += ctx.lepage*w.staple(dxdw[mu], mu, nu);
        )
        HISQLOOP2(
          dsi = w.staple(dxdw[nu], nu, i);
          HISQLOOP3(
            sj = ctx.c3*w.staple(nu, j);
            dsnu += ctx.c2*dsi + ctx.c3*w.staple(dsi, nu, j);
            cnu += w.staple(ctx.c2*dxdw[mu] + ctx.c3*w.staple(dxdw[mu], mu, j), mu, i);
          )
          snu += w.staple(ctx.c2*w.link(nu) + sj, nu, i);
          dxdu[mu] += w.stapleDerivative(sj, dsi, mu, nu);
        )
        dxdu[mu] += w.stapleDerivative(snu, dxdw[nu], mu, nu); // Eqn 3a
        dxdu[mu] += w.stapleDerivative(dsnu, mu, nu);          // Eqn 3a
        dxdu[mu] += w.staple(cnu, mu, nu);                     // Eqn 3b
    ) )

    HISQLOOP0(dxdu[mu] = adj(dxdu[mu])) // correct for orientation of staples
    HISQNAIK(
      GaugeLorentzField dxdwww = toLorentz(cell.Exchange(dXdWWW));
      HISQLOOP0(
        si = w.Cshift(w.link(mu), mu, FORWARD);
        sj = si*w.Cshift(si, mu, FORWARD)*adj(dxdwww[mu]);
        snu = sj;
        for (int term = 0; term < 2; ++term) {
          sj = w.Cshift(adj(si)*sj*w.link(mu), mu, BACKWARD);
          snu += sj;
        }
        dxdu[mu] += ctx.naik*snu;
    ) )
    dXdU = w.toTightGrid(toGauge(dxdu));
  }

  void smearDerivative(
    GaugeField& dXdU,
    const GaugeField& dXdW,
    const GaugeField& dXdWWW,
    const GaugeField& W
  ) { 
    HISFContext asqtadCtx(ASQL1, ASQL3, ASQL5, ASQL7, LEPAGE, NAIK); 
    smearDerivative(dXdU, dXdW, dXdWWW, W, asqtadCtx); 
  }

  void smearDerivative(
    GaugeField& dXdU,
    const GaugeField& dXdW,
    const GaugeField& W,
    const HISFContext ctx
  ) { smearDerivative(dXdU, dXdW, dXdW, W, ctx); }

  void smearDerivative(GaugeField& dXdU, const GaugeField& dXdW, const GaugeField& W) { 
    HISFContext fat7Ctx(F7L1, F7L3, F7L5, F7L7);  
    smearDerivative(dXdU, dXdW, W, fat7Ctx);
  }

  void projectionDerivative(
    GaugeField& dVdU, 
    const GaugeField& dZdV, 
    const GaugeField& U,
    const HISFContext ctx
  ) {
    UnitaryProjection<Gimpl> projection(
      ctx.eigenCutoff, 
      ctx.backupSVD, 
      ctx.svdTolerance
    );
    projection.derivative(dVdU, dZdV, U); 
  }

  void projectionDerivative(GaugeField& dVdU, const GaugeField& dZdV, const GaugeField& U) {
    UnitaryProjection<Gimpl> projection(
      REUNITDERIVCUTOFF, 
      BACKUPSVD, 
      BACKUPSVDTOLERANCE
    );
    projection.derivative(dVdU, dZdV, U);
  }

  void projectionDerivative(
    GaugeField& dVdU, 
    const GaugeField& dZdV, 
    const GaugeField& V,
    const GaugeField& U
  ) { // projection constructor inputs don't matter because no projection performed
    UnitaryProjection<Gimpl> projection(1e-20, true, 1e-8);
    projection.derivative(dVdU, dZdV, V, U); 
  }

};

// undefine macros to prevent conflicts
#undef HISQLOOP0
#undef HISQLOOP1
#undef HISQLOOP2
#undef HISQLOOP3

NAMESPACE_END(Grid);

#endif // QCD_UTILS_HISQ_IMPL_H