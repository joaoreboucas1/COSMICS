# `COSMICS`: Cosmological Initial Conditions and Microwave Anisotropy Codes

> [!Note]
> This repository contains the `COSMICS` code, written by Edmund Bertschinger and Paul Bode. I do not claim any ownership of the code in this repository. All rights are reserved to the original authors, see `LICENSE`.

`COSMICS` is one of the first Einstein-Boltzmann codes, written by Edmund Bertschinger and Paul Bode. It implements the Einstein-Boltzmann system of linear equations described in the seminal paper [Cosmological Perturbation Equations in the Synchronous and Conformal Newtonian Gauges](https://arxiv.org/pdf/astro-ph/9506072). The documentation can be found at [this Arxiv paper](https://arxiv.org/pdf/astro-ph/9506070). However, the website that used to host this code is now offline. The code can still be found at [Ed Bertschinger's website](https://web.mit.edu/edbert/). I decided to put the code in Github for historical reference.

## Building

While the repository contains the original code, to build the programs in a modern Linux system, the user should modify the first two options in `Make_files/Make.LINUX` as follows:
```
#
# Compiler options for Linux.
#
F77 = gfortran            # Original is g77
F77FLAGS = -std=legacy -O # Original is just -O
...
```

Afterwards, just compile the code normally with the command `make linux`.

## References

- https://arxiv.org/pdf/astro-ph/9506070
- https://arxiv.org/pdf/astro-ph/9506072
- https://web.mit.edu/edbert/
- Search on Wayback Machine for http://arcturus.mit.edu/cosmics/