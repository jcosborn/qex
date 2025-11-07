# QEX: Quantum EXpressions lattice field theory framework

## Description of QEX fork:

This is my personal fork of the Quantum EXpressions lattice field theory framework. I use it both for production running of project that I am currently involved in and for testing new ideas out.

- The code under [src/stagg_pv_hmc](https://github.com/ctpeterson/qex/tree/devel/src/stagg_pv_hmc) was used for a determination of the Nf = 8 and 12 beta-function. It contains both software for Hamiltonian Monte Carlo and gradient flow. It is outdated and I don't recommend using it for production running, as I no longer support it. Nonetheless, it should run efficiently; see the instructions below for a description of how to get it running. 

- The code under [src/mcmc](https://github.com/ctpeterson/qex/tree/devel/src/mcmc) was my attempt to build a comprehensive, yet simple to use, software suite for running Hamiltonian Monte Carlo simulations of staggered fermions. See [src/examples/fxpvyhmc.nim](https://github.com/ctpeterson/qex/blob/devel/src/examples/fxpvyhmc.nim) for an example of how it can be used. It deploys both unsmeared and nHYP-smeared staggered fermions with or without rooting; however, the unsmeared parts of the code are not tested. If you'd like to use this for production, please contact me via email so that I can help you check to make sure it is running that way it should for your specific purpose, particularly with respect to the integrators, which must be run with specific constraints to ensure their reliability. My contact information is below. The multi-shift solver was written by me and can be found at [src/solvers/cgm.nim](https://github.com/ctpeterson/qex/blob/devel/src/examples/fxpvyhmc.nim). I'd like to return to refactoring this code one day.

- Software for unrooted highly-improved staggered quark simulations can be found under [src/alphas](https://github.com/ctpeterson/qex/blob/devel/src/examples/fxpvyhmc.nim). It also contains a nice gradient flow suite.

- I've various tests of Nambu Hamiltonian Monte Carlo under [src/nhmc](https://github.com/ctpeterson/qex/blob/devel/src/examples/fxpvyhmc.nim).

- There are a number of less interesting/useful bits of code scattered throughout. This is the power of the design of QEX: it makes testing out wild ideas incredibly simple.

Contact: curtistaylorpetersonwork@gmail.com

## Old description of QEX fork:
Modular code for fast deployment of Hybrid Monte Carlo (HMC) with many options for simulating with/without nHYP smeared (or unsmeared) staggered fermions and Pauli-Villars (PV) bosons. Most options (gauge action, boundary conditions, nHYP smearing parameters, and much more) are provided in an XML file, such that the action being simulated can be quickly/easily modified without needing to modify the source code and/or recompile. Source code can be found [here](https://github.com/ctpeterson/qex/tree/devel/src/stagg_pv_hmc). See [input_hmc.xml](https://github.com/ctpeterson/qex_staghmc/blob/devel/src/stagg_pv_hmc/input_hmc.xml) for a complete list of tunable parameters for the gauge-fermion-PV action, HMC specifications, etc. 

For deployment, make sure that you 1.) follow the installation instructions given in the base QEX respository and 2.) install [MDevolve](https://github.com/jxy/MDevolve) using Nimble. Once you've done this, simply run "make staghmc_spv" in your build directory and you're done. If you're confused because the whole process seemed far too easy, I assure you that it really was that easy and encourage you to thank James Osborn and Xiao-Yong Jin for the time & effort that they've put into QEX. Instructions for running the code can be found at the top of [staghmc_spv.nim](https://github.com/ctpeterson/qex_staghmc/blob/devel/src/stagg_pv_hmc/staghmc_spv.nim).

This fork of QEX also ships with a modular gauge flow code, allowing for Wilson, adjoint-plaquette, and any variation of rectangular action for the flow. The gauge flow code can be found [here](https://github.com/ctpeterson/qex/tree/devel/src/flow).

There is also an XY model code with Wolff cluster updates that can be found [here](https://github.com/ctpeterson/qex/tree/devel/src/xy_cluster_mc). The XY model simulation code also comes with options for performing U(1) gradient flow with the XY model action as the flow action, which is equivalent to a gradient flow of the XY model action with an explicit constraint that preserves the unit norm of the spins on the lattice. 

For the U(1) gradient flow of the XY model action, I use [Arraymancer](https://mratsim.github.io/Arraymancer/index.html) for vectorization. The latter gradient flow code is not as efficient as the gauge flow code that is native to QEX, and it would be desireable to have a version of the U(1) gradient flow code that is build within the QEX framework. One option is simply to compile the [gauge flow code](https://github.com/ctpeterson/qex/tree/devel/src/flow) in this repository with --nc=1 and to modify the action to be the XY model action; this may come in a future version of this repository, but it is unlikely.

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
