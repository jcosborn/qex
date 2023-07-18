# QEX: Quantum EXpressions lattice field theory framework

## Description of QEX fork:
Modular code for fast deployment of Hybrid Monte Carlo (HMC) with many options for simulating with/without nHYP smeared (or unsmeared) staggered fermions and Pauli-Villars (PV) bosons. Most options (gauge action, boundary conditions, nHYP smearing parameters, and much more) are provided in an XML file, such that the action being simulated can be quickly/easily modified without needing to modify the source code and/or recompile. Source code can be found [here](https://github.com/ctpeterson/qex/tree/devel/src/stagg_pv_hmc). See input_hmc.xml for a complete list of tunable parameters for the gauge-fermion-PV action, HMC specifications, and much more. 

For deployment, make sure that you 1.) follow the installation instructions given in the base QEX respository and 2.) install [MDevolve](https://github.com/jxy/MDevolve) using Nimble. Once you've done this, simply run "make staghmc_spv" in your build directory and you're done. If you're confused because the whole process seemed far too easy, I assure you that it really was that easy and encourage you to thank James Osborn and Xiao-Yong Jin for the time & effort that they've put into QEX. 

This fork of QEX also ships with a modular gauge flow code, allowing for Wilson, Adjoint-Plaquette, and any variation of rectangular action for the flow. The gauge flow code can be found [here](https://github.com/ctpeterson/qex/tree/devel/src/flow).

There is also an XY model code with Wolff cluster updates that can be found [here](https://github.com/ctpeterson/qex/tree/devel/src/xy_cluster_mc). The XY model simulation code also comes with options for performing U(1) gradient flow with the XY model action as the flow action, which is equivalent to a gradient flow of the XY model action with an explicit constraint that preserves the unit norm of the spins on the lattice. 

For the U(1) gradient flow of the XY model action, I use [Arraymancer](https://mratsim.github.io/Arraymancer/index.html) for vectorization. The latter gradient flow code is not as efficient as the gauge flow code that is native to QEX, and it would be desireable to have a version of the U(1) gradient flow code that is build within the QEX framework.

## Description of QEX:
QEX is a high-level framework for lattice field operations
written in the language [Nim](https://nim-lang.org).

It provides optimized lattice field operations, including SIMD support,
for CPU architectures (native GPU support is currently experimental).
Since Nim compiles to native C/C++, directly calling any C/C++ lattice
code or library from QEX is relatively easy to do.

Some simple code examples are here
 [ex0.nim](src/examples/ex0.nim)
 [ex1.nim](src/examples/ex1.nim).

It currently supports
- U(1), SU(2..4) gauge fields in any dimension
- SciDAC I/O
- Gauge fixing
- Staggered solver and forces (Asqtad, HISQ, nHYP)
- Wilson solver (no clover yet)
- Interface for Chroma, Grid, QUDA interoperability

Installation guide: [INSTALL.md](INSTALL.md)

Build guide: [BUILD.md](BUILD.md)

Further examples:
- [tests/examples](tests/examples)
