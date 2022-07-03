# pmns_phi64

This repository contains standard C source code for the PMNS systems with $\phi=2^{64}$.

In order to perform the tests, it is necessary to follow the this procedure, in a shell:

- execute the `./rdpmc.sh` script
- go into the directory of the desired version
- type `make -B`
- then `./main`

The Makefile specifies `gcc-10`, which is the gcc version we used in our tests. Previous versions may not be able to compile the source code as it is.
