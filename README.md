# QTL

This is a collection of (small) test data, algorithms and functions about 
breeding and genetics.  They are from many of my repos, mainly on GitHub.
These datawere used many times, and hence majorly confirmd.

They are grouped in the following sections:
- `Data`, majorly for test purpose.
- `Simulation`, to simulate various scenarios in breeding
- `Matrices`, matrix manipulation algorithms
- `Notes`, including some standard procedures
- to be continued


# Version log
Newest first

- v0.2.0: Genome scan on 30M SNP.
  From this version, an `Xps` file is only included in `src/xps.jl` in the 
  corresponding version. This is to avoid function incompatible as the package develops.
  Later, a minor version increment means that a new test was done. The test was done
  in a branch, which was later merged into main, resulting an increment of minor version.
- v0.1.2: (large) GRM construction and inverse
  - QTL.MAT.grm
- v0.1.1: quick genotype simulation, no LD.
  - QTL.MISC.quick_g, matrix I/O
