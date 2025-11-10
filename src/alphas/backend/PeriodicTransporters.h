/**
 * @file PeriodicTranspoerters.h
 * @brief Interface for periodic gauge transporters
 * @author Curtis Taylor Peterson
 */

#pragma once 

#ifndef QCD_UTILS_PERIODIC_TRANSPORTERS_H
#define QCD_UTILS_PERIODIC_TRANSPORTERS_H

#include <Grid/Grid.h>

NAMESPACE_BEGIN(Grid);

//
// macros
//

// make scope more explicit
#define ACCELERATOR_SCOPE(exec) { exec; } \

// shorten call to get stencil entry in declaration
#define NEW_STENCIL_ENTRY(se, st, mu, n)              \
  GeneralStencilEntry const *se = st.GetEntry(mu, n); \

// shorten call to set stencil entry
#define SET_STENCIL_ENTRY(se, st, mu, n) se = st.GetEntry(mu, n); \

// shorten call to coalesced read in function
#define ACCREAD(u, se)                                         \
  coalescedReadGeneralPermute(u[se->_offset], se->_permute, Nd) \

// shorten call to coalesced write
#define ACCWRITE(wu, u) coalescedWrite(wu, u) \

//
// convenient data structures
//

enum TransportHeading {FORWARD = 0, BACKWARD = Nd};

//
// periodic transporter
//

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
  int depth;
  int mu;
  std::unique_ptr<GaugeLinkField> _ubuf;
  std::unique_ptr<GaugeLinkField> _vbuf;
  std::shared_ptr<GeneralLocalStencil> _sbuf;

public:
  PeriodicTransporter(){}

  PeriodicTransporter(
    std::shared_ptr<GeneralLocalStencil> stencil, 
    const GaugeLinkField& U, 
    int mu,
    int depth
  ): _sbuf(stencil), mu(mu), depth(depth) {
    _vbuf = std::make_unique<GaugeLinkField>(GaugeLinkField(U.Grid()));
    _ubuf = std::make_unique<GaugeLinkField>(U);
  }

public:
  /** @brief application of gauge transporter to operand field */
  inline const GaugeLinkField CovShift(
    const GaugeLinkField& v, 
    TransportHeading heading
  ) {
    ACCELERATOR_SCOPE(
      GeneralLocalStencilView sbuf_v = (*_sbuf).View(AcceleratorRead);
      autoView(v_v, v, AcceleratorRead);
      autoView(ubuf_v, (*_ubuf), AcceleratorRead);
      autoView(vbuf_v, (*_vbuf), AcceleratorWrite);
      auto forward = [&](int n) {
        NEW_STENCIL_ENTRY(se, sbuf_v, mu, n);
        ACCWRITE(vbuf_v[n], ubuf_v[n]*ACCREAD(v_v, se));
      };
      auto backward = [&](int n) {
        NEW_STENCIL_ENTRY(se, sbuf_v, mu + BACKWARD, n);
        ACCWRITE(vbuf_v[n], adj(ACCREAD(ubuf_v, se))*ACCREAD(v_v, se));
      };
      if (heading == FORWARD)
      { accelerator_for(n, v_v.size(), Simd::Nsimd(), forward(n);); }
      else accelerator_for(n, v_v.size(), Simd::Nsimd(), backward(n););
    )
    return (*_vbuf);
  }

  /** @brief application of shift operation to operand field */
  inline const GaugeLinkField Cshift(const GaugeLinkField& v, TransportHeading heading) {
    ACCELERATOR_SCOPE(
      GeneralLocalStencilView sbuf_v = (*_sbuf).View(AcceleratorRead);
      autoView(v_v, v, AcceleratorRead);
      autoView(vbuf_v, (*_vbuf), AcceleratorWrite);
      accelerator_for(n, v_v.size(), Simd::Nsimd(), {
        NEW_STENCIL_ENTRY(se, sbuf_v, mu + heading, n); 
        ACCWRITE(vbuf_v[n], ACCREAD(v_v, se));
      });
    )
    return (*_vbuf);
  }

  /** @brief link elongation of link buffer */
  inline const GaugeLinkField elongate() {
    ACCELERATOR_SCOPE(
      GeneralLocalStencilView sbuf_v = (*_sbuf).View(AcceleratorRead);
      autoView(ubuf_v, (*_ubuf), AcceleratorRead);
      autoView(vbuf_v, (*_vbuf), AcceleratorWrite);
      accelerator_for(n, ubuf_v.size(), Simd::Nsimd(), {
        auto u = ubuf_v[n];
        for (int d = 0; d < depth; ++d) {
          NEW_STENCIL_ENTRY(se, sbuf_v, mu + d*Nd, n); 
          u = u*ACCREAD(ubuf_v, se);
        }
        ACCWRITE(vbuf_v[n], u);
      });
    )
    return (*_vbuf);
  }

  /** @brief derivative of link buffer elongation */
  inline const GaugeLinkField elongationDerivative(const GaugeLinkField& dsdwww) {
    assert(depth == 2 && "elongation derivative only implemented for Naik");
    ACCELERATOR_SCOPE(
      GeneralLocalStencilView sbuf_v = (*_sbuf).View(AcceleratorRead);
      autoView(dsdwww_v, dsdwww, AcceleratorRead);
      autoView(ubuf_v, (*_ubuf), AcceleratorRead);
      autoView(vbuf_v, (*_vbuf), AcceleratorWrite);
      accelerator_for(n, ubuf_v.size(), Simd::Nsimd(), {
        NEW_STENCIL_ENTRY(se1p, sbuf_v, mu, n);             // +1
        NEW_STENCIL_ENTRY(se1m, sbuf_v, mu + 2*Nd, n);      // -1
        NEW_STENCIL_ENTRY(se2p, sbuf_v, mu + Nd, n);        // +2
        NEW_STENCIL_ENTRY(se2m, sbuf_v, mu + Nd + 2*Nd, n); // -2
        ACCWRITE(
          vbuf_v[n],
          ACCREAD(ubuf_v, se1p)*ACCREAD(ubuf_v, se2p)*dsdwww_v[n] + \
          ACCREAD(ubuf_v, se1p)*ACCREAD(dsdwww_v, se1m)*ACCREAD(ubuf_v, se1m) + \
          ACCREAD(dsdwww_v, se2m)*ACCREAD(ubuf_v, se2m)*ACCREAD(ubuf_v, se1m)
        );
      });
    )
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
  inline void set_vbuf(const GaugeLinkField& v) 
  { _vbuf.reset(); _vbuf = std::make_unique<GaugeLinkField>(v); }

  /** @brief for modifying "w" buffer */
  inline void set_wbuf(const GaugeLinkField& w) 
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
  std::unique_ptr<PaddedCell> _pcell;
  std::unique_ptr<GaugeLinkField> _vbuf;
  std::shared_ptr<GeneralLocalStencil> _sbuf;
  Transporter _t[Nd];

private:
  GaugeLinkField toLink(const GaugeField& U, int mu)
  { return PeekIndex<LorentzIndex>(U, mu); }

public:
  PeriodicTransporters(PaddedCell& pcell, const GaugeField& Uin) {
    _pcell = std::make_unique<PaddedCell>(pcell);

    int depth = _pcell->depth;
    auto U = _pcell->ExchangePeriodic(Uin);
    auto* grid = U.Grid();
    std::vector<Coordinate> shifts(2*depth*Nd, 0);

    for (int mu = 0; mu < Nd; ++mu) { 
      for (int d = 0; d < depth; ++d) {
        shifts[mu + d*Nd][mu] = d + 1; 
        shifts[mu + d*Nd + depth*Nd][mu] = -(d + 1); 
    } }

    _sbuf = std::make_shared<GeneralLocalStencil>(GeneralLocalStencil(grid, shifts));
    _vbuf = std::make_unique<GaugeLinkField>(GaugeLinkField(grid));

    for (int mu = 0; mu < Nd; ++mu) 
      _t[mu] = Transporter(_sbuf, toLink(U, mu), mu, depth);
  }


public:
  /** @brief wrapped halo exchange */
  inline const GaugeField toPaddedGrid(const GaugeField& U) 
  { return _pcell->ExchangePeriodic(U); }

  /** @brief wrapped extraction from padding */
  inline const GaugeField toTightGrid(const GaugeField& U) 
  { return _pcell->Extract(U); }

  /** @brief wrapped extraction from padding */
  inline const GaugeLinkField toTightGrid(const GaugeLinkField& U) 
  { return _pcell->Extract(U); }

public:
  /** @brief cartesian shift (only periodic) */
  inline GaugeLinkField Cshift(const GaugeLinkField& u, int mu, TransportHeading heading) 
  { return _t[mu].Cshift(u, heading); }

public:
  /** @brief calculate symmetric staple: without mu pre-shift */
  inline const GaugeLinkField staple(
    const GaugeLinkField& v, // mu link
    const GaugeLinkField& u, // nu link
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
    auto* grid = v.Grid();
    GaugeLinkField us(grid), ls(grid);
    
    us = Zero();
    ls = Zero();
    ACCELERATOR_SCOPE( // scope: calculate staple
      GeneralLocalStencilView sbuf_v = (*_sbuf).View(AcceleratorRead);
      ACCELERATOR_SCOPE( // calculate upper staple and unshifted lower staple
        autoView(u_v, u, AcceleratorRead);
        autoView(v_v, v, AcceleratorRead);
        autoView(us_v, us, AcceleratorWrite);
        autoView(ls_v, ls, AcceleratorWrite);
        accelerator_for(n, us_v.size(), Simd::Nsimd(), {
          NEW_STENCIL_ENTRY(se_mu, sbuf_v, mu, n);
          NEW_STENCIL_ENTRY(se_nu, sbuf_v, nu, n);
          auto su_v = ACCREAD(u_v, se_mu);
          ACCWRITE(us_v[n], u_v[n]*ACCREAD(v_v, se_nu)*adj(su_v));
          ACCWRITE(ls_v[n], adj(u_v[n])*v_v[n]*su_v);
        });
      )
      ACCELERATOR_SCOPE( // add upper staple to downward-shifted lower staple 
        autoView(us_v, us, AcceleratorRead);
        autoView(ls_v, ls, AcceleratorRead);
        autoView(vbuf_v, (*_vbuf), AcceleratorWrite);
        accelerator_for(n, us_v.size(), Simd::Nsimd(), {
          NEW_STENCIL_ENTRY(se, sbuf_v, nu + BACKWARD, n);
          ACCWRITE(vbuf_v[n], us_v[n] + ACCREAD(ls_v, se));
        });
      ) 
    )
    return (*_vbuf);
  }

  /** @brief calculate symmetric staple: without mu pre-shift */
  inline const GaugeLinkField staple(const GaugeLinkField& v, int mu, int nu) 
  { return staple(v, link(nu), mu, nu); }

  /** @brief calculate symmetric staple using buffer fields: without mu pre-shift **/
  inline const GaugeLinkField staple(int mu, int nu) 
  { return staple(link(mu), link(nu), mu, nu); }

public:
  inline const GaugeLinkField stapleDerivative(
    const GaugeLinkField& v, // mu link
    const GaugeLinkField& u, // nu link
    const GaugeLinkField& c, // chain
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
    auto* grid = v.Grid();
    GaugeLinkField uds(grid), lds(grid);

    uds = Zero();
    lds = Zero();
    ACCELERATOR_SCOPE( // calculate staple derivative
      GeneralLocalStencilView sbuf_v = (*_sbuf).View(AcceleratorRead);
      ACCELERATOR_SCOPE( // upper contribution and unshifted lower contribution
        autoView(v_v, v, AcceleratorRead);
        autoView(u_v, u, AcceleratorRead);
        autoView(c_v, c, AcceleratorRead);
        autoView(lds_v, lds, AcceleratorWrite);
        autoView(uds_v, uds, AcceleratorWrite);
        accelerator_for(n, lds_v.size(), Simd::Nsimd(), {
          NEW_STENCIL_ENTRY(se_mu, sbuf_v, mu, n);
          NEW_STENCIL_ENTRY(se_nu, sbuf_v, nu, n);
          auto su_v = ACCREAD(u_v, se_mu);
          auto sc_v = ACCREAD(c_v, se_mu);
          auto sv_v = ACCREAD(v_v, se_nu);
          ACCWRITE(uds_v[n], c_v[n]*sv_v*adj(su_v) + u_v[n]*sv_v*adj(sc_v));
          ACCWRITE(lds_v[n], adj(c_v[n])*v_v[n]*su_v + adj(u_v[n])*v_v[n]*sc_v);
        });
      )
      ACCELERATOR_SCOPE( // shift lower contribution and to result
        autoView(lds_v, lds, AcceleratorRead);
        autoView(uds_v, uds, AcceleratorRead);
        autoView(vbuf_v, (*_vbuf), AcceleratorWrite);
        accelerator_for(n, lds_v.size(), Simd::Nsimd(), {
          NEW_STENCIL_ENTRY(se, sbuf_v, nu + BACKWARD, n);
          ACCWRITE(vbuf_v[n], uds_v[n] + ACCREAD(lds_v, se));
        });
      )
    )
    return (*_vbuf);
  }

  /** @brief nu-oriented symmetric staple derivative w/o passing of middle link */
  inline const GaugeLinkField stapleDerivative(
    const GaugeLinkField& v,
    const GaugeLinkField& c,
    int mu,
    int nu
  ) { return stapleDerivative(link(mu), v, c, mu, nu); }

  /** @brief nu-oriented symmetric staple derivative w/o explicit middle/side links */
  inline const GaugeLinkField stapleDerivative(const GaugeLinkField& c, int mu, int nu) 
  { return stapleDerivative(link(mu), link(nu), c, mu, nu); }

public:
  /** @brief calculation of long link */
  inline const GaugeLinkField elongate(int mu, int length) {
    assert(Nd == 4 && "link elongation only implemented in 4D");
    return _t[mu].CovShift(_t[mu].CovShift(link(mu), FORWARD), FORWARD);
  }

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

// undefine macros to prevent conflicts
#undef ACCELERATOR_SCOPE
#undef NEW_STENCIL_ENTRY
#undef SET_STENCIL_ENTRY
#undef ACCREAD
#undef ACCWRITE

NAMESPACE_END(Grid);

#endif // QCD_UTILS_PERIODIC_TRANSPORTERS_H