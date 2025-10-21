## QEX: Quantum EXpressions lattice field theory framework

QEX is a high-level framework for lattice field operations
written in the language [Nim](https://nim-lang.org).

It provides optimized lattice field operations, including SIMD support,
for CPU architectures (native GPU support is currently experimental).
Since Nim compiles to native C/C++, directly calling any C/C++ lattice
code or library from QEX is relatively easy to do.

Some simple code examples are here
 [ex0.nim](src/examples/ex0.nim)
 [ex1.nim](src/examples/ex1.nim).

Some of the operations it supports are
- U(1), SU(2..4) gauge fields in any dimension
- SciDAC I/O
- Gauge fixing
- Staggered solver and forces (Asqtad, HISQ, nHYP)
- Wilson solver (no clover yet)
- Interface for Chroma, Grid, QUDA interoperability

QEX is under active development and new capabilities are constantly being added.

<!--SKIP
All content before this sign is replicated in the Introduction chapter of nimibook documentation
-->

Installation guide: [INSTALL.md](INSTALL.md)

Build guide: [BUILD.md](BUILD.md)

Further examples:
- [tests/examples](tests/examples)
