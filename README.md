# RBR
---
RBR (Row-By-Row) is a fast solver for community detection type problems. It is mainly
implemented in C language, with OpenMP support. We also provide more friendly
interfaces for MATLAB users.

## Introduction
RBR is a modularity-based community detection method, which solves the 
modularity maximization problem

![mm](_img/mm.jpg)

via non-convex relaxation:

![non-convex](_img/non-convex.jpg)

The matrix `U` is optimized row by row asynchronously.

RBR requires the number of clusters `k` to perform clustering. It has a good
performance even if `k` is small.

## How to Install
RBR has only been tested on Linux currently. To use RBR you will need a recent
Linux operating system with:
- A recent C compiler with OpenMP support (gcc >= 4.8)
- A C++ compiler that supports C++11 (g++ >= 4.8 is fine)
- GNU Make
- CBLAS or Intel MKL

First download the source code and unzip anywhere you like.

Then change into the directory `src` and edit `Makefile` with your favorite
text editor. Basically, you will need to specify the compilers and the location
of the required libraries.

### Compile with Intel MKL
It is highly recommended to use Intel MKL as external BLAS routines,
which is publicly available at [here](https://software.intel.com/en-us/mkl).

To compile with Intel MKL, one should set `MKL` variable in `Makefile` to
`on`, and specify `MKLROOT` (this can be done using Intel MKL's initialization
script). 

### Compile with CBLAS
If no Intel MKL is available, RBR can be linked against any implementation of
CBLAS (the performance will be lower, though). To use CBLAS, one should set
`MKL` variable to `off`, and specify the name of the CBLAS header and libraries.
For example:
```
CBLAS_HEADER = cblas.h
CBLAS_LIBS = -lblas
```

### Compilation
After you have done with the configuration part of Intel MKL/CBLAS, simply
type `make` in the directory `src`. If no errors occur, you should see
the output files:
- `libDCRBR.a`: static library of RBR
- `rbr`: ELF executable of RBR

### Compile the MATLAB interface
To compile the MATLAB interface you'll need a recent MATLAB distribution.
First add MATLAB binary directory to `PATH`. Then type
```
make mex
```
to compile MEX files. If there are no errors, you should see the output
file:
- `matlab/mex_rbr.mexa64`: MEX file (64-bit)

## Usage
To use RBR, type `./rbr -h` to see the usage. You can also run the examples
we provide:
```
./rbr -v --full ../examples/polblogs 2
./rbr -v ../examples/amazon 100 5
```

RBR currently only supports [Rutherford Boeing (RB) Sparse Matrix File Format]
(http://people.math.sc.edu/Burkardt/data/rb/rb.html)
as input. More formats will be supported in the future release.

## References
- [Junyu Zhang, Haoyang Liu, Zaiwen Wen, and Shuzhong Zhang, A sparse completely positive relaxation of the modularity maximization for community detection](https://epubs.siam.org/doi/abs/10.1137/17M1141904)

## The authors
We hope that the package is useful for your applications. Please feel free
to contact the authors if you have any comments or bug reports.

- Haoyang Liu (liuhaoyang@pku.edu.cn)
- Zaiwen Wen  (wenzw@pku.edu.cn)

## Copyright
RBR

Copyright (C) 2019  Haoyang Liu (liuhaoyang@pku.edu.cn)
                    Zaiwen Wen  (wenzw@pku.edu.cn)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.