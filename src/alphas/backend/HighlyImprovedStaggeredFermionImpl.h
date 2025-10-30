/*************************************************************************************
Grid physics library, www.github.com/paboyle/Grid

Source file: ./lib/qcd/utils/HighlyImprovedStaggeredFermionImpl.h

Copyright (C) 2015

Author: Curtis Taylor Peterson <curtistaylorpetersonwork@gmail.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

/**
 * @file HighlyImprovedStaggeredFermionImpl.h
 * @brief Interface for implementation of highly improved staggered fermions (HISQ)
 * @author Curtis Taylor Peterson and David Clarke
 * @details
 * This header file is meant to act as an interface for both Grid and the 
 * MILC codebase to utilize the "highly improved staggered quark" action
 * within Grid. The interface itself is composed of three core component
 * classes:
 * 
 * * PeriodicTransporter: Optimized gauge transport with PaddedCell
 * 
 * * PeriodicTransporters: Container of GaugeTransporter objects for 
 *   multi-directional gauge transport
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
#include "UnitaryProjection.h"

NAMESPACE_BEGIN(Grid);

//
// useful consts
//

// usual 4D HISQ coefficients: defined for convenience
const int // Lepage & Naik
  LEPAGE = -1.0/8.0,
  NAIK = -1.0/24.0;
const int // fat-7
  F7L1 = -LEPAGE,
  F7L3 = 0.25*F7L1,
  F7L5 = 0.25*F7L3,
  F7L7 = -0.125*NAIK;
const int // asqtad
  ASQL1 = -8.0*LEPAGE,
  ASQL3 = F7L3,
  ASQL5 = F7L5,
  ASQL7 = F7L7; 

//
// macros
//

// shorten call to get stencil entry in declaration
#define NEW_STENCIL_ENTRY(se, st, mu, n)              \
  GeneralStencilEntry const *se = st.GetEntry(mu, n); \

// shorten call to set stencil entry
#define SET_STENCIL_ENTRY(se, st, mu, n) se = st.GetEntry(mu, n); \

// shorten call to coalesced read in function
#define HISQREAD(u, se)                                         \
  coalescedReadGeneralPermute(u[se->_offset], se->_permute, Nd) \

// shorten call to coalesced write
#define HISQWRITE(wu, u) coalescedWrite(wu, u) \

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
#define HISQLEPAGE(exec) if (lepage != 0.0) { exec; } \

// encapsulates conditional execution of Naik calculations for cleanliness
#define HISQNAIK(exec) if (naik != 0.0) { exec; } \

//
// staggered implementation parameters
//

struct StagImplParams {
  Coordinate dirichlet; // Blocksize of dirichlet BCs
  int  partialDirichlet;
  AcceleratorVector<Complex,Nd> boundary_phases;
  StagImplParams()
  {partialDirichlet = 0; dirichlet.resize(0); boundary_phases.resize(Nd,1.0);};
  StagImplParams(std::vector<Complex> phi): boundary_phases(phi)
  {partialDirichlet = 0; dirichlet.resize(0);};
};

//
// periodic transporter
//

enum HEADING {FORWARD = 0, BACKWARD = Nd};

template <class Gimpl>
class PeriodicTransporter: public Gimpl
{
/**
 * @brief Object representing *periodic* gauge transporter
 * @author Curtis Taylor Peterson
 * @details
 * This is just about the most naive way that one could implement a class that
 * uses the padded cell to represent a gauge transporter. Many optimizations are
 * possible, so take this as the minimal example that must be improved upon. The
 * advantage of this over regular covariant shifts in Grid (for periodic gauge
 * fields) is this defines operations that do nothing more than read/write to a
 * buffer. Not having to create temporaries lends to a massive speedup.
 */
public: INHERIT_GIMPL_TYPES(Gimpl)

private:
  int mu;
  std::unique_ptr<GaugeLinkField> _ubuf;
  std::unique_ptr<GaugeLinkField> _vbuf;
  std::shared_ptr<GeneralLocalStencil> _sbuf;

public:
  PeriodicTransporter(){}

  PeriodicTransporter(
    std::shared_ptr<GeneralLocalStencil> stencil, 
    const GaugeLinkField &U, 
    int mu
  ): mu(mu), _sbuf(stencil) {
    _vbuf = std::make_unique<GaugeLinkField>(GaugeLinkField(U.Grid()));
    _ubuf = std::make_unique<GaugeLinkField>(U);
  }

public:
  /** @brief application of gauge transporter to operand field */
  inline const GaugeLinkField CovShift(const GaugeLinkField &v, HEADING heading) {
    {
      GeneralLocalStencilView sbuf_v = (*_sbuf).View(AcceleratorRead);
      autoView(v_v, v, AcceleratorRead);
      autoView(ubuf_v, (*_ubuf), AcceleratorRead);
      autoView(vbuf_v, (*_vbuf), AcceleratorWrite);
      auto forward = [&](int n) {
        NEW_STENCIL_ENTRY(se, sbuf_v, mu, n);
        HISQWRITE(vbuf_v[n], ubuf_v[n]*HISQREAD(v_v, se));
      };
      auto backward = [&](int n) {
        NEW_STENCIL_ENTRY(se, sbuf_v, mu + BACKWARD, n);
        HISQWRITE(vbuf_v[n], adj(HISQREAD(v_v, se))*HISQREAD(ubuf_v, se));
      };
      if (heading == FORWARD)
      { accelerator_for(n, v_v.size(), Simd::Nsimd(), forward(n);); }
      else accelerator_for(n, v_v.size(), Simd::Nsimd(), backward(n););
    }
    return (*_vbuf);
  }

  /** @brief application of shift operation to operand field */
  inline const GaugeLinkField Cshift(const GaugeLinkField &v, HEADING heading) {
    {
      GeneralLocalStencilView sbuf_v = (*_sbuf).View(AcceleratorRead);
      autoView(v_v, v, AcceleratorRead);
      autoView(vbuf_v, (*_vbuf), AcceleratorWrite);
      accelerator_for(n, v_v.size(), Simd::Nsimd(), {
        NEW_STENCIL_ENTRY(se, sbuf_v, mu + heading, n); 
        HISQWRITE(vbuf_v[n], HISQREAD(v_v, se));
      });
    }
    return (*_vbuf);
  }

public:
  /** @brief return copy of input buffer */
  inline const GaugeLinkField link() {return (*_ubuf);}

  /** @brief return copy of output buffer */
  inline const GaugeLinkField field() {return (*_vbuf);}

  /** @brief return copy of general local stencil */
  inline const GeneralLocalStencil stencil() {return (*_sbuf);}

public:
  /** @brief for modifying "v" buffer */
  inline void set_vbuf(const GaugeLinkField &v) 
  { _vbuf.reset(); _vbuf = std::make_unique<GaugeLinkField>(v); }

  /** @brief for modifying "w" buffer */
  inline void set_wbuf(const GaugeLinkField &w) 
  { _ubuf.reset(); _ubuf = std::make_unique<GaugeLinkField>(w); }
};

//
// periodic transporter container
//

template <class Gimpl>
class PeriodicTransporters: public Gimpl
{
/**
 * @brief PeriodicTransporter container
 * @author Curtis Taylor Peterson
 * @details
 * Naive container for saving and indexing a collection of PeriodicTransporters. 
 * Note of caution: all periodic gauge transporters use the same buffer for the 
 * stencil and the outputs.
 */

public: INHERIT_GIMPL_TYPES(Gimpl)

private:
  typedef PeriodicTransporter<Gimpl> Transporter;

private:
  std::unique_ptr<GaugeLinkField> _vbuf;
  std::shared_ptr<GeneralLocalStencil> _sbuf;
  Transporter _t[Nd];

private:
  GaugeLinkField _get(const GaugeField &U, int mu)
  { return PeekIndex<LorentzIndex>(U, mu); }

public:
  PeriodicTransporters(PaddedCell &pcell, const GaugeField &Uin) {
    auto U = pcell.Exchange(Uin);
    auto *grid = U.Grid();
    std::vector<Coordinate> shifts(2*Nd, 0);

    for (int mu = 0; mu < Nd; ++mu){shifts[mu][mu] = 1; shifts[mu + Nd][mu] = -1;}
    _sbuf = std::make_shared<GeneralLocalStencil>(GeneralLocalStencil(grid, shifts));
    _vbuf = std::make_unique<GaugeLinkField>(GaugeLinkField(grid));

    for (int mu = 0; mu < Nd; ++mu){_t[mu] = Transporter(_sbuf, _get(U, mu), mu);}
  }

public:
  /** @brief cartesian shift (only periodic) */
  inline GaugeLinkField Cshift(const GaugeLinkField &u, int mu, HEADING heading) 
  { return _t[mu].Cshift(u, heading); }

public:
  /** @brief calculate symmetric staple: without mu pre-shift */
  inline const GaugeLinkField staple(
    const GaugeLinkField &v, // mu link
    const GaugeLinkField &u, // nu link
    int mu, 
    int nu
  ) {
    /**
     * @brief Calculates "symmetric" staple (i.e., staple + its reflection)
     * @author Curtis Taylor Peterson
     * @details
     * Calculates closed/symmetric staple with orientation
     *         ---🠢
     *         🠡   🠣
     * ν       x---x
     * 🠡       🠣   🠡
     * -🠢 μ    ---🠢 
     * where leftmost "x" is at "n" and rightmost "x" is
     * at n + μ. As such "v = Umu" and ""
     */ 
    auto *grid = v.Grid();
    GaugeLinkField us(grid), ls(grid);
    
    us = Zero();
    ls = Zero();
    { // scope: calculate staple
      GeneralLocalStencilView sbuf_v = (*_sbuf).View(AcceleratorRead);

      { // scope: calculate upper staple and unshifted lower staple
        autoView(u_v, u, AcceleratorRead);
        autoView(v_v, v, AcceleratorRead);
        autoView(us_v, us, AcceleratorWrite);
        autoView(ls_v, ls, AcceleratorWrite);
        accelerator_for(n, us_v.size(), Simd::Nsimd(), {
          NEW_STENCIL_ENTRY(se_mu, sbuf_v, mu, n);
          NEW_STENCIL_ENTRY(se_nu, sbuf_v, nu, n);
          auto su_v = HISQREAD(u_v, se_mu);
          HISQWRITE(us_v[n], u_v[n]*HISQREAD(v_v, se_nu)*adj(su_v));
          HISQWRITE(ls_v[n], adj(u_v[n])*v_v[n]*su_v);
        });
      }

      { // scope: add upper staple to downward-shifted lower staple 
        autoView(us_v, us, AcceleratorRead);
        autoView(ls_v, ls, AcceleratorRead);
        autoView(vbuf_v, (*_vbuf), AcceleratorWrite);
        accelerator_for(n, us_v.size(), Simd::Nsimd(), {
          NEW_STENCIL_ENTRY(se, sbuf_v, nu + BACKWARD, n);
          HISQWRITE(vbuf_v[n], us_v[n] + HISQREAD(ls_v, se));
        });
    } 
  }
    return (*_vbuf);
  }

  /** @brief calculate symmetric staple: without mu pre-shift */
  inline const GaugeLinkField staple(const GaugeLinkField &v, int mu, int nu) 
  { return staple(v, link(nu), mu, nu); }

  /** @brief calculate symmetric staple using buffer fields: without mu pre-shift **/
  inline const GaugeLinkField staple(int mu, int nu) 
  { return staple(link(mu), link(nu), mu, nu); }

  /** @brief calculate symmetric staple: with mu pre-shift **/
  inline const GaugeLinkField staple(
    const GaugeLinkField &v, 
    const GaugeLinkField &psu, 
    int nu
  ) {
    /**
     * @brief Calculates symmetrci staple with u_nu preshifted in mu-direction
     * @author Curtis Taylor Peterson
     * @details See symmetric staple method above for details
     */
    {
      GeneralLocalStencilView sbuf_v = (*_sbuf).View(AcceleratorRead);
      autoView(v_v, v, AcceleratorRead);
      autoView(u_v, link(nu), AcceleratorRead);
      autoView(psu_v, psu, AcceleratorRead);
      autoView(vbuf_v, (*_vbuf), AcceleratorWrite);
      accelerator_for(n, v_v.size(), Simd::Nsimd(), {
        NEW_STENCIL_ENTRY(se, sbuf_v, nu, n);
        auto s = u_v[n]*HISQREAD(v_v, se)*adj(psu_v[n]); 
        SET_STENCIL_ENTRY(se, sbuf_v, nu + BACKWARD, n);
        HISQWRITE(
          vbuf_v[n], 
          s + adj(HISQREAD(u_v, se))*HISQREAD(v_v, se)*HISQREAD(psu_v, se)
        );
      });
    }
    return (*_vbuf);
  }

  /** @brief calculate symmetric staple using buffer fields: with mu pre-shift **/
  inline const GaugeLinkField staple(const GaugeLinkField &v, int nu) 
  { return staple(v, field(nu), nu); }

public:
  inline const GaugeLinkField stapleDerivative(
    const GaugeLinkField &v, // mu link
    const GaugeLinkField &u, // nu link
    const GaugeLinkField &c, // chain
    int mu,
    int nu
  ) {
    /**
     * @brief Calculates nu-ordiented derivative of symmetric staple
     * @author Curtis Taylor Peterson
     * @details
     * Calculates nu-oriented derivative of closed/symmetric staple. See diagram in
     * symmetric staple method for orientation of symmetric staple. Derivative "D"
     * of just nu links (what this method calculates) is
     *            ----🠢      ----🠢 
     *            ⮾    🠣     🠡   ⮾ 
     * ν     D =  x----x  +  x----x
     * 🠡          ⮾    🠡     🠣   ⮾
     * -🠢 μ       ----🠢      ----🠢 
     *             [1a]       [1b] 
     * where replacement of "🠣" or "🠡" with ⮾ indicates replacement of link with 
     * contribution from chain rule (i.e., action of derivative).
     */
    auto *grid = v.Grid();
    GaugeLinkField uds(grid), lds(grid);

    uds = Zero();
    lds = Zero();
    { // scope: calculate staple derivative
      GeneralLocalStencilView sbuf_v = (*_sbuf).View(AcceleratorRead);

      { // scope: upper contribution and unshifted lower contribution
        autoView(v_v, v, AcceleratorRead);
        autoView(u_v, u, AcceleratorRead);
        autoView(c_v, c, AcceleratorRead);
        autoView(lds_v, lds, AcceleratorWrite);
        autoView(uds_v, uds, AcceleratorWrite);
        accelerator_for(n, lds_v.size(), Simd::Nsimd(), {
          NEW_STENCIL_ENTRY(se_mu, sbuf_v, mu, n);
          NEW_STENCIL_ENTRY(se_nu, sbuf_v, nu, n);
          auto su_v = HISQREAD(u_v, se_mu);
          auto sc_v = HISQREAD(c_v, se_mu);
          auto sv_v = HISQREAD(v_v, se_nu);
          HISQWRITE(uds_v[n], c_v[n]*sv_v*adj(su_v) + u_v[n]*sv_v*adj(sc_v));
          HISQWRITE(lds_v[n], adj(c_v[n])*v_v[n]*su_v + adj(u_v[n])*v_v[n]*sc_v);
        });
      }

      { // scope: shift lower contribution and to result
        autoView(lds_v, lds, AcceleratorRead);
        autoView(uds_v, uds, AcceleratorRead);
        autoView(vbuf_v, (*_vbuf), AcceleratorWrite);
        accelerator_for(n, lds_v.size(), Simd::Nsimd(), {
          NEW_STENCIL_ENTRY(se, sbuf_v, nu + BACKWARD, n);
          HISQWRITE(vbuf_v[n], uds_v[n] + HISQREAD(lds_v, se));
        });
      }
    }
    return (*_vbuf);
  }

  /** @brief nu-oriented symmetric staple derivative w/o passing of middle link */
  inline const GaugeLinkField stapleDerivative(
    const GaugeLinkField &v,
    const GaugeLinkField &c,
    int mu,
    int nu
  ) { return stapleDerivative(link(mu), v, c, mu, nu); }

  /** @brief nu-oriented symmetric staple derivative w/o explicit middle/side links */
  inline const GaugeLinkField stapleDerivative(const GaugeLinkField &c, int mu, int nu) 
  { return stapleDerivative(link(mu), link(nu), c, mu, nu); }

public:
  /** @brief index transporters -- mutable */
  inline Transporter &operator[](int mu) {return _t[mu];}

  /** @brief index transporters -- not mutable */
  inline const Transporter &operator[](int mu) const {return _t[mu];}

public:
  /** @brief return copy of input buffer */
  inline const GaugeLinkField link(int mu) {return _t[mu].link();}

  /** @brief return copy of output buffer */
  inline const GaugeLinkField field(int mu) {return _t[mu].field();}

  /** @brief return copy of general local stencil */
  inline const GeneralLocalStencil stencil(int mu) {return _t[mu].stencil();}
};

//
// highly improved staggered fermion implementation
//

template <class Gimpl>
class HighlyImprovedStaggeredFermionImpl: Gimpl {
  
public: INHERIT_GIMPL_TYPES(Gimpl);

private:
  typedef typename Simd::scalar_type GridScalar;
  typedef typename std::vector<GaugeLinkField> GaugeLorentzField;
  typedef typename std::vector<ComplexField> StaggeredPhases;

public:
  GridBase *grid;
  PaddedCell cell;

  StagImplParams params;
  StaggeredPhases stagPhases;

  bool backupSVD;
  RealD projectionCutoff;
  RealD svdTolerance;

public:
  HighlyImprovedStaggeredFermionImpl(
    GridCartesian *grid, 
    const StagImplParams &params,
    RealD projectionCutoff = 1e-8,
    bool backupSVD = true,
    RealD svdTolerance = 5e-5
  ):
    cell(1, grid), 
    stagPhases(Nd, grid), 
    params(params), 
    backupSVD(backupSVD),
    projectionCutoff(projectionCutoff),
    svdTolerance(svdTolerance) {
    assert(Nc == 3 && "HISQ only suppored for SU(3)");
    assert(Nd == 4 && "HISQ only supported for 4 dimensions");
    this->grid = (GridBase*)cell.grids.back();
    getStaggeredPhases(grid, stagPhases);
  }

  HighlyImprovedStaggeredFermionImpl(GridCartesian *grid):
    cell(1, grid), 
    stagPhases(Nd, grid), 
    params(StagImplParams({1, 1, 1, -1})), 
    backupSVD(true),
    projectionCutoff(1e-8),
    svdTolerance(5e-5) {
    assert(Nc == 3 && "HISQ only suppored for SU(3)");
    assert(Nd == 4 && "HISQ only supported for 4 dimensions");
    this->grid = (GridBase*)cell.grids.back();
    getStaggeredPhases(grid, stagPhases);
  }

private:
  void getStaggeredPhases(GridBase *grid, StaggeredPhases &Phases) {
    /**
     * @brief HISQ gauge configuration constructor
     * @author Curtis Taylor Peterson, Peter Boyle
     * @details
     * Staggered phases follow the "MILC convention", which treats the fourth 
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
    Lattice<iScalar<vInteger>> x(grid), y(grid), t(grid);
    Lattice<iScalar<vInteger>> tx(grid), txy(grid), xyzt(grid), coor(grid);
    ComplexField phases(grid);

    LatticeCoordinate(x, 0);
    LatticeCoordinate(y, 1);
    LatticeCoordinate(t, 3);
    tx = t + x;
    txy = tx + y;

    for (int mu = 0; mu < Nd; mu++) { 
      // for boundary phases
      int N = grid->GlobalDimensions()[mu] - 1;
      auto bpha = params.boundary_phases[mu];
      GridScalar dirichlet(real(bpha), imag(bpha));

      // staggered phases x boundary phases
      LatticeCoordinate(coor, mu);
      phases = 1.0;
      if (mu == 0){phases = where(mod(t, 2) == (Integer)0,   phases, -phases);}
      if (mu == 1){phases = where(mod(tx, 2) == (Integer)0,  phases, -phases);}
      if (mu == 2){phases = where(mod(txy, 2) == (Integer)0, phases, -phases);}
      Phases[mu] = where(coor == (Integer)N, dirichlet*phases, phases);
    }
  }

private:
  // wrap PokeIndex (insert gauge link field into gauge field layout)
  inline void insertLink(GaugeField &Uout, const GaugeLinkField &U, int mu)
  { PokeIndex<LorentzIndex>(Uout, U, mu); }

  // wrap PeekIndex (extract gauge link field form gauge field layout)
  inline const GaugeLinkField toLink(const GaugeField &U, int mu)
  { return PeekIndex<LorentzIndex>(U, mu); }

  // break up memory layout of gauge field into std::vector of gauge link fields
  inline const GaugeLorentzField toLorentz(const GaugeField &U) 
  { GaugeLorentzField u(Nd, grid); HISQLOOP0(u[mu] = toLink(U, mu);) return u; }

  // aggregate std::vector of gauge link fields into layout of gauge field
  inline const GaugeField toGauge(GaugeLorentzField &u) 
  { GaugeField U(grid); HISQLOOP0(insertLink(U, u[mu], mu);) return U; }

public:
  /** @brief rephases gauge links with staggered and Dirichlet boundary phases */
  void rephase(GaugeField &X, const GaugeField &W) { 
    for (int mu = 0; mu < Nd; ++mu) 
      insertLink(X, stagPhases[mu]*toLink(W, mu), mu); 
  }

public:
  void smear(
    GaugeField &X,
    const GaugeField &W, 
    std::vector<RealD> coeffs,
    RealD lepage = 0.0
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
     * Then an Asqtad link Vμ(n) in four dimensions can be expressed as a 
     * composition of sums
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

    RealD c0 = coeffs[0], c1 = coeffs[1], c2 = coeffs[2], c3 = coeffs[3];
    PeriodicTransporters<Gimpl> w(cell, W);
    GaugeLorentzField x(Nd, grid);
    GaugeLinkField s3(grid), s5(grid); 
    
    HISQLOOP0(                                             //
      x[mu] = (c0 - 2.0*(double)(Nd-1)*lepage)*w.link(mu); // Eqn [5a] + lepage fix
      HISQLOOP1(                                           //
        s3 = w.staple(w.link(mu), mu, nu);                 // 
        x[mu] += c1*s3;                                    // Eqn [5b]
        HISQLEPAGE(x[mu] += lepage*w.staple(s3, mu, nu))   // lepage
        HISQLOOP2(                                         // 
          s5 = w.staple(s3, mu, i);                        // 
          x[mu] += c2*s5;                                  // Eqn [5c]
          HISQLOOP3(x[mu] += c3*w.staple(s5, mu, j))       // Eqn [5d]
    ) ) )                                                  //
    X = cell.Extract(toGauge(x));
  }

  void project(GaugeField &W, const GaugeField &U) { 
    UnitaryProjection<Gimpl> projection(0.0, backupSVD, svdTolerance);
    projection.project(W, U); 
  } 

public:
  void smearDerivative(
    GaugeField &dXdU,
    const GaugeField &dXdW,
    const GaugeField &dXdWWW,
    const GaugeField &W,
    std::vector<RealD> coeffs,
    RealD lepage = 0.0,
    RealD naik = 0.0
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
     * I'm going here, we need to consider these contributions "from the perspective
     * of the force". In other words, before we take the derivative, we think of the 
     * the contribution from the chain rule ∂S/∂Vδ as completing a plaquette 
     * (just with the like ∂S/∂Vδ pointing in the "wrong" direction); in Fig. 1, this 
     * is the dotted line below the staple. We then think of the contribution from 
     * the derivative on the left and top link by "rotating" the perspective in 
     * Fig. 1, such that the link that the derivative ∂/∂Uμ(n) hits lies on the 
     * x-axis where the contribution from ∂S/∂Vδ was. We then get
     *                 Uμ(n)
     * ν      🠠----   ---🠢 
     * 🠡      ⮾   🠡 + 🠡   🠣 (Fig 2)
     * -🠢 μ   ---🠢    -⮾-
     *        Uμ(n)
     *        [3a]     [3b]
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
     * from the 3-link staples in a similar manner. I'll close off by noting that 
     * this way of visualizing the Wirtinger derivative works generally; in fact,
     * it's an incredibly useful tool for quickly reasoning about the structure of
     * any Hamiltonian (hybrid) Monte Carlo force. 
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

    // TODO:
    // - use periodic exchange for both main HISQ loop & Naik???

    RealD c0 = coeffs[0], c1 = coeffs[1], c2 = coeffs[2], c3 = coeffs[3];
    PeriodicTransporters<Gimpl> w(cell, W);
    GaugeLorentzField dxdw = toLorentz(cell.Exchange(dXdW));
    GaugeLorentzField dxdu(Nd, grid);
    GaugeLinkField cnu(grid), ci(grid);
    GaugeLinkField snu(grid), si(grid), sj(grid);
    GaugeLinkField dsnu(grid), dsi(grid);

    // fat7 + lepage (lepage won't execute if lepage == 0.0)
    HISQLOOP0(
      dxdu[mu] = (c0 - 6.0*lepage)*dxdw[mu];
      HISQLOOP1(
        snu = c1*w.link(nu);
        cnu = c1*dxdw[mu];
        HISQLEPAGE(
          snu += lepage*w.staple(w.link(nu), nu, mu);
          cnu += lepage*w.staple(dxdw[mu], mu, nu);
          dsnu = lepage*w.staple(dxdw[nu], nu, mu);
        )
        if (lepage == 0.0) dsnu = Zero();
        HISQLOOP2(
          dsi = w.staple(dxdw[nu], nu, i);
          HISQLOOP3(
            sj = c3*w.staple(nu, j);
            dsnu += c2*dsi + c3*w.staple(dsi, nu, j);
            cnu += w.staple(c2*dxdw[mu] + c3*w.staple(dxdw[mu], mu, j), mu, i);
          )
          snu += w.staple(c2*w.link(nu) + sj, nu, i);
          dxdu[mu] += w.stapleDerivative(sj, dsi, mu, nu);
        )
        dxdu[mu] += w.stapleDerivative(snu, dxdw[nu], mu, nu);
        dxdu[mu] += w.stapleDerivative(dsnu, mu, nu);
        dxdu[mu] += w.staple(cnu, mu, nu);
      )
      dxdu[mu] = adj(dxdu[mu]);
    )

    // naik (won't execute if naik == 0.0)
    // move this outside of this method?
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
        dxdu[mu] += naik*snu;
    ) )

    // extract from padded layout and return
    dXdU = adj(cell.Extract(toGauge(dxdu)));
  }

  void smearDerivative(
    GaugeField &dXdU,
    const GaugeField &dXdW,
    const GaugeField &W,
    std::vector<RealD> coeffs,
    RealD lepage = 0.0
  ) { smearDerivative(dXdU, dXdW, dXdW, W, coeffs, lepage, 0.0); }

  void projectionDerivative(
    GaugeField &dVdU, 
    const GaugeField &dZdV, 
    const GaugeField &U
  ) { 
    UnitaryProjection<Gimpl> projection(projectionCutoff, backupSVD, svdTolerance);
    projection.derivative(dVdU, dZdV, U); 
  }

public:
  /** @brief constructs long link (Naik) */
  void elongate(GaugeField &WWW, const GaugeField &W, RealD naik = 1.0) { 
    PeriodicTransporters<Gimpl> w(cell, W);
    GaugeLorentzField www(Nd, grid);
    HISQLOOP0(www[mu] = w[mu].CovShift(w[mu].CovShift(w.link(mu), FORWARD), FORWARD);)
    WWW = cell.Extract(toGauge(www));
    if (naik != 0.0) WWW *= naik;
  }

};

// undefine macros
#undef NEW_STENCIL_ENTRY
#undef SET_STENCIL_ENTRY
#undef HISQREAD
#undef HISQWRITE
#undef HISQLOOP0
#undef HISQLOOP1
#undef HISQLOOP2
#undef HISQLOOP3

NAMESPACE_END(Grid);

#endif // QCD_UTILS_HISQ_IMPL_H