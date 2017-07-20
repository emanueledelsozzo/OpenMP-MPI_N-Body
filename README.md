# OpenMP-MPI_N-Body

Author: 
Emanuele Del Sozzo (emanuele.delsozzo@polimi.it)

This repository contains a parallelization of N-Body simulation algorithm using OpenMP and MPI libraries.

This code was developed for the "Parallel Programming With MPI and OpenMP" Ph.D. class held at Politecnico di Milano in collaboration with Cineca (https://www.cineca.it/en).

The repository is organized as follows:
- "mpi_code" folder contains 3 different parallelizations of N-Body simulation using MPI
- "original_code" folder contains the original code of N-Body simulation provided by Cineca
- "openmp_code" folder contains a parallelization of N-Body simulation using OpenMP
- "presentation" folder contains the .pptx file of the project presentation
- "profiling" folder contains the profiling of the application performed using Valgrind toolchain

The N-Body simulation application was parallelized on the "Galileo" cluster (http://www.hpc.cineca.it/hardware/galileo).