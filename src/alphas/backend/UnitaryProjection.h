/**
 * @file UnitaryProjection.h
 * @brief Defines object representing unitary projection
 * @author Curtis Taylor Peterson
 * @details
 * This header file defines an object for performing unitary projection of gauge
 * links. The projection is performed using Cayley-Hamilton as a default and 
 * optionally Jacobi-based singular value decomposition (SVD) as a fallback.
 * 
 * References:
 * * Hasenfratz, A. et al.: https://doi.org/10.1103/PhysRevD.78.014515
 * * MILC Collaboration: https://doi.org/10.1103/PhysRevD.75.054502
 * * Quantum EXpressions (QEX): https://github.com/jcosborn/qex
 * * QOPQDP [SciDAC]: https://github.com/usqcd-software/qopqdp
 * 
 * Acknowledgements:
 *   Curtis Taylor Peterson would like to thank James Osborn and Xiao-Yong Jin for 
 *   developing/testing the implementations of unitary projection in Quantum 
 *   EXpressions, from which this implementation is based and has been tested against.
 *   
 *   This material is based upon work supported by the U.S. Department of Energy, 
 *   Office of Science, Office of Advanced Scientific Computing Research, Scientific 
 *   Discovery through Advanced Computing (SciDAC) program.
*/

#pragma once

#ifndef QCD_UTILS_UNITARYPROJECTION_H
#define QCD_UTILS_UNITARYPROJECTION_H

#include <Grid/Grid.h>

NAMESPACE_BEGIN(Grid);

const RealD SMALL = std::numeric_limits<double>::epsilon();

template<class Gimpl> 
class UnitaryProjection: public Gimpl {

public: INHERIT_GIMPL_TYPES(Gimpl);

private:
  typedef typename Simd::scalar_type GridScalar;
  typedef iScalar<iScalar<iMatrix<GridScalar, Nc>>> GridScalarMatrix;
  
  typedef typename Eigen::Matrix<ComplexD, Nc, Nc> EigenScalarMatrix;
  typedef typename Eigen::JacobiSVD<EigenScalarMatrix> EigenSVD;

  RealD cutoff, svdtol;
  bool backupSVD;

public:
  UnitaryProjection(
    RealD cutoff = SMALL, 
    bool backupSVD = false,
    RealD svdtol = 1e-8
  ): cutoff(cutoff), svdtol(svdtol), backupSVD(backupSVD) 
  { assert(Nc == 3 && "unitary projection only supported for Nc = 3 for now"); }

private:
  EigenScalarMatrix toEigen(const GridScalarMatrix& u) {
    EigenScalarMatrix eu;
    for (int i = 0; i < Nc; ++i) {
      for (int j = 0; j < Nc; ++j) {
        GridScalar uij = u()()(i, j);
        eu(i, j) = ComplexD(real(uij), imag(uij));
    } }
    return eu;
  }

  GridScalarMatrix toGrid(const EigenScalarMatrix& u) {
    GridScalarMatrix gu;
    for (int i = 0; i < Nc; ++i) {
      for (int j = 0; j < Nc; ++j) {
        ComplexD uij = u(i, j);
        gu()()(i, j) = GridScalar(real(uij), imag(uij));
    } }
    return gu;
  }

private:
  void _adjugate3(GaugeLinkField &Ai, const GaugeLinkField& A) {
    GridBase *grid = A.Grid();
    GaugeLinkField T(grid);
    LatticeComplex trA(grid),trA2(grid);
    T = A*A;
    trA = trace(A), trA2 = trace(T);
    Ai = T - trA*A;
    T = 1.0;
    Ai += 0.5*(trA*trA - trA2)*T;
  }

  void _inverse3(GaugeLinkField &Ai, const GaugeLinkField& A)
  { _adjugate3(Ai,A); Ai = Ai/Determinant(A); }

  void _sylvester3(
    GaugeLinkField& X, 
    const GaugeLinkField& A, 
    const GaugeLinkField& C
  ) {
    GridBase *grid = A.Grid();
    GaugeLinkField adjA(grid);
    GaugeLinkField AC(grid), CA(grid), ACA(grid);
    GaugeLinkField adjAC(grid), CadjA(grid), adjACadjA(grid);
    LatticeComplex unit(grid), t(grid), s(grid), r(grid);
    LatticeComplex c0(grid), c1(grid), c2(grid), c3(grid);
    _adjugate3(adjA,A);
    t = trace(A), s = trace(adjA);
    r = peekColour(A,0,0)*peekColour(adjA,0,0);
    r += peekColour(A,0,1)*peekColour(adjA,1,0);
    r += peekColour(A,0,2)*peekColour(adjA,2,0);
    AC = A*C;
    CA = C*A;
    ACA = AC*A;
    adjAC = adjA*C;
    CadjA = C*adjA;
    adjACadjA = adjAC*adjA;
    unit = 1.0;
    c2 = 0.5*unit/(s*t - r);
    c0 = c2*(s + t*t);
    c3 = c2*t;
    c1 = c3/r;
    X = c0*C + c1*adjACadjA; 
    X += c2*(ACA - adjAC - CadjA);
    X -= c3*(AC + CA);
  }

  LatticeComplex _absmin(const LatticeComplex& x, const LatticeComplex& y) { 
    LatticeReal xr = toReal(x);
    LatticeReal yr = toReal(y);
    xr = abs(xr);
    return where(xr <= yr, x, y); 
  }

  LatticeComplex _absmax(const LatticeComplex& x, const LatticeComplex& y) { 
    LatticeReal xr = toReal(x);
    LatticeReal yr = toReal(y);
    xr = abs(xr);
    return where(xr >= yr, x, y); 
  }

  void _bound(LatticeComplex& x) {
    LatticeReal xr = toReal(x);
    x = where(xr < cutoff, x + cutoff, x);
  }

  void _eigs3(
    LatticeComplex& f0,
    LatticeComplex& f1,
    LatticeComplex& f2,
    const GaugeLinkField& q,
    const GaugeLinkField& q2
  ) {
    GridBase *grid = q.Grid();
    Complex k1 = 1.0/3.0, k2 = 0.5*k1, k3 = 2.0*M_PI*k1;
    LatticeComplex ir(grid), uv(grid);
    LatticeComplex a0(grid), a1(grid), a2(grid);

    ir = SMALL, uv = 1.0; 

    a0 = k1*real(trace(q));
    a1 = k2*real(trace(q2)); 
    a2 = k2*real(trace(q*q2));

    f1 = a0*a0;
    a2 += a0*(f1 - 3.0*a1);
    a1 -= 0.5*f1;
    a1 = sqrt(abs(a1));

    a2 = _absmin(a2/_absmax(a1*a1*a1, ir), uv);
    a1 *= 2.0;
    a2 = k1*acos(a2);

    f0 = a0;
    f1 = f0 + a1*cos(a2);
    f2 = f0 + a1*cos(a2 + k3);
    f0 += a1*cos(a2 - k3);
    _bound(f0), _bound(f1), _bound(f2);
  }

private:
  void _projectU3(GaugeLinkField& v, const GaugeLinkField& u) {
    /**
     * @brief U(3) unitary projection via Cayley-Hamilton or SVD
     * @author Curtis Taylor Peterson
     * @details
     * This method implements a U(3) projection of a general complex 3x3 matrix 
     * using Cayley-Hamilton or the Jacobi singular value decomposition implemented
     * by Eigen. For details about the Cayley-Hamilton approach, see the references
     * provided above; namely the OG paper by Hasenfratz et al and later work by the 
     * MILC collaboration. Please note that this method is modelled after the approach
     * taken in Quantum EXpressions by James Osborn and Xiao-Yong Jin.
     */
    GridBase *grid = u.Grid();
    GaugeLinkField unity(grid), q(grid), q2(grid);
    LatticeComplex e0(grid), e1(grid), e2(grid);
    LatticeComplex f0(grid), f1(grid), f2(grid);
    LatticeComplex unit(grid), detA(grid), detB(grid);

    // Cayley-Hamilton: eigenvalues of q = u†u
    unit = 1.0, unity = 1.0;
    q = adj(u)*u; 
    q2 = q*q;
    _eigs3(e0, e1, e2, q, q2);
    detA = Determinant(q), detB = e0*e1*e2;

    // Cayley-Hamilton: "u, v, w" coefficients [Eqn. C6 of PRD(75)054502]
    f0 = sqrt(e0), f1 = sqrt(e1), f2 = sqrt(e2);
    e0 = f0 + f1 + f2;
    e1 = f0*f1;
    e2 = e1*f2;
    e1 += f0*f2 + f1*f2;

    // Cayley-Hamilton: "f0, f1, f2" coefficients [Eqn. C7 of PRD(75)054502]
    f2 = e2*(e0*e1 - e2);
    f2 = unit/f2;
    f1 = e0*e0;
    f0 = e0*e1*e1 - e2*(f1 + e1);
    f0 *= f2;
    f1 = e0*(2.0*e1 - f1) - e2;
    f1 *= f2;
    f2 *= e0;
    v = u*(f0*unity + f1*q + f2*q2);

    // Jacobi-based singular value decomposition: fallback for ill-conditioned links
    // conditions for falling back on SVD: https://doi.org/10.1103/PhysRevD.75.054502
    if (backupSVD) {{
      autoView(detA_v, detA, CpuRead);
      autoView(detB_v, detB, CpuRead);
      autoView(u_v, u, CpuRead);
      autoView(e0_v, e0, CpuRead);
      autoView(e1_v, e1, CpuRead);
      autoView(e2_v, e2, CpuRead);
      autoView(v_v, v, CpuWrite);
      thread_for(n, grid->lSites(), { // TODO: mask
        bool detDiffTooLarge, e0TooSmall, e1TooSmall, e2TooSmall;
        Coordinate lcoor;
        GridScalar localDetA, localDetB;
        GridScalar locale0, locale1, locale2;

        grid->LocalIndexToLocalCoor(n, lcoor);
        peekLocalSite(localDetA, detA_v, lcoor);
        peekLocalSite(localDetB, detB_v, lcoor);
        peekLocalSite(locale0, e0_v, lcoor);
        peekLocalSite(locale1, e1_v, lcoor);
        peekLocalSite(locale2, e2_v, lcoor);

        detDiffTooLarge = abs(localDetA - localDetB) > svdtol;
        e0TooSmall = abs(locale0) < svdtol;
        e1TooSmall = abs(locale1) < svdtol;
        e2TooSmall = abs(locale2) < svdtol;

        if (detDiffTooLarge or e0TooSmall or e1TooSmall or e2TooSmall) {
          GridScalarMatrix gu;
          EigenScalarMatrix eu, ev = EigenScalarMatrix::Zero();
          
          peekLocalSite(gu, u_v, lcoor);
          EigenSVD svd(toEigen(gu), Eigen::ComputeFullU | Eigen::ComputeFullV);
          ev = svd.matrixU() * svd.matrixV().adjoint();
          pokeLocalSite(toGrid(ev), v_v, lcoor);
        }
      });
    }}
  }

  void _derivativeU3(
    GaugeLinkField& dvdu, 
    const GaugeLinkField& dzdv,
    const GaugeLinkField& v,
    const GaugeLinkField& u
  ) {
    /**
     * @brief Derivative of unitary projection
     * @author Curtis Taylor Peterson
     * @details
     * This is a very clever approach to calculating the derivative of the unitary 
     * projection; it was  invented by the authors of the original HISQ paper 
     * [PRD72(2007)054502] and expanded upon by James Osborn and Xiao-Yong Jin in the 
     * Quantum EXpressions code. The objective is to calculate 
     * (1) dQ/dU = CZ + U dZ/dU,
     * where
     * (2) Z = (X'X)^{-1/2}
     * (3) Q = XZ,
     * (4) C = dX/dU.
     * We can solve for dZ/dU using the identity
     * (5) dZ/dU = -Q'CZ = Y dZ/dU + dZ/dU Y
     * with
     * (6) Y = (X'X)^{1/2}.
     * This is a so-called "Sylvester" system of equations, and it has an analytic 
     * solution for N = 3 that is calculated in the _sylvester3 method.
     */
    GridBase *grid = u.Grid();
    GaugeLinkField t1(grid), t2(grid), t3(grid);

    _inverse3(t1, v);        // (u'u)^(1/2) u^-1
    t2 = t1*u;               // (u'u)^(1/2)  [6]
    _inverse3(t3, t2);       // (u'u)^(-1/2) [3]
    dvdu = dzdv*t3;          //
    t1 = adj(v)*dvdu;        // second equality of Eqn. [5]
    _sylvester3(t3, t2, t1); // solve Sylvester [5]
    t2 = t3 + adj(t3);       // d/dX & d/dX'
    dvdu -= u*t2;            //
  }

public:
  void project(GaugeLinkField& v, const GaugeLinkField& u) { _projectU3(v, u); }

  void project(GaugeField& V, const GaugeField& U) {
    GaugeLinkField v(U.Grid());
    for (int mu = 0; mu < Nd; ++mu) {
      _projectU3(v, PeekIndex<LorentzIndex>(U, mu));
      PokeIndex<LorentzIndex>(V, v, mu);
    }
  }

  void derivative(GaugeField& dVdU, const GaugeField& dZdV, const GaugeField& U) {
    GaugeLinkField dvdu(U.Grid());
    GaugeLinkField u(U.Grid()), v(U.Grid());
    for (int mu = 0; mu < Nd; ++mu){
      u = PeekIndex<LorentzIndex>(U, mu);
      _projectU3(v, u);
      _derivativeU3(dvdu, PeekIndex<LorentzIndex>(dZdV, mu), v, u);
      PokeIndex<LorentzIndex>(dVdU, dvdu, mu);
    }
  }

  void derivative(
    GaugeField& dVdU, 
    const GaugeField& dZdV, 
    const GaugeField& V, 
    const GaugeField& U
  ) {
    // be careful w/ using this method: you want to make sure 
    // that cutoff used here is consistent w/ cutoff used to get v
    GaugeLinkField dvdu(U.Grid());
    for (int mu = 0; mu < Nd; ++mu){
      _derivativeU3(
        dvdu, 
        PeekIndex<LorentzIndex>(dZdV, mu), 
        PeekIndex<LorentzIndex>(V, mu),
        PeekIndex<LorentzIndex>(U, mu)
      );
      PokeIndex<LorentzIndex>(dVdU, dvdu, mu);
    }
  }

};

NAMESPACE_END(Grid);

#endif // QCD_UTILS_UNITARYPROJECTION_H