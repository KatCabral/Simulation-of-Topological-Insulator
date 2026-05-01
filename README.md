# Simulation-of-Topological-Insulator
Author: Katarina Cabral

**Introduction**
-------------------
The study of topological materials is a relatively new and interesting branch of condensed matter physics. The unique properties of topologically non-trivial devices have great potential in the field of quantum informatics and for the development of robust, energy-efficient electronic and spintronic devices. This is the software for a project which investigated the mechanisms governing topological phase transitions induced by external perturbations of applied pressure for a three-dimensional cross-correlated topological material in which the spin-splitting Rashba effect occurs. The software undertakes numerical tight-binding calculations, simulating the topological phase transition in the giant Rashba semiconductor BiTeI.

**Code Function**
-------------------
This program, Wannier_band_structure, computes the band structure evolution across a simulated topological phase transition in the giant Rashba semiconductor BiTeI. It interpolates between "trivial" and "topological" Hamiltonians to visualise how band inversion and surface states evolve under adiabatic pressure application.

The code models a topological phase transition induced by hydrostatic pressure in BiTeI. Using tight-binding Hamiltonians generated from maximally localised Wannier functions, it linearly interpolates between trivial and topological data to create an intermediate “band touching point” (BTP) phase.

The Hamiltonian is constructed in k-space via Fourier transform from real-space hopping parameters contained in BiTeI_hr_triv.dat and BiTeI_hr_top.dat. 

**Data File Strucuture**
--------------------------
**BiTeI_hr_triv.dat and BiTeI_hr_top.dat**

The first line contains a single integer giving the total number of Wannier functions in the model. For example, a value of 18 means that each Hamiltonian matrix `H(R)` is 18 by 18.

The second line gives the total number of distinct lattice translation vectors R included in the model. For example, 1155 indicates that 1155 different R-vectors are stored in the file.

The following section lists one integer per R-vector, giving the degeneracy (number of equivalent vectors related by symmetry) for each lattice translation. These integers are read sequentially, even though they may be written in rows for compactness.

After the degeneracy list, the file contains the Hamiltonian matrix elements. Each line corresponds to one complex Hamiltonian element `H_ij(R)` and has the form:

        Rx Ry Rz i j Re(Hij) Im(Hij)

`Rx`, `Ry`, and `Rz` are the integer components of the lattice vector `R`. `i` and `j` are the indices of the Wannier functions, ranging from 1 up to the total number given at the start of the file (18 in our example). `Re(Hij)` and `Im(Hij)` are the real and imaginary parts of the Hamiltonian matrix element connecting orbital `i` in the home cell to orbital j in the cell displaced by `R`. These `H_ij(R)` values represent the hopping parameters of the system: the strength and phase of the coupling between Wannier orbitals `i` and `j` separated by lattice vector `R`.

The Hamiltonian is Hermitian, which means `H_ij(R) = H_ji*(-R)`.

Each line therefore represents one element of the Hamiltonian for a particular lattice vector.

**BiTeI.nnkp**

Real-space lattice vectors: Between begin real_lattice and end real_lattice, three lines list the primitive vectors of the crystal in angstroms, each with x, y, z components.

Reciprocal lattice vectors: Between begin recip_lattice and end recip_lattice, three lines list the reciprocal lattice vectors in inverse angstroms, used to convert real-space positions to k-space.

K-points: Between begin kpoints and end, the first line gives the total number of k-points. Each following line lists one k-point as fractional coordinates `(kx, ky, kz)` in the reciprocal lattice. These points define where the Hamiltonian `H(k)` will be evaluated.

The file is read sequentially by the Fortran program to extract lattice geometry, k-points, and calculation options for constructing `H(k)` from the real-space Hamiltonian `H(R)`.

**Code Descripton and Running**
---------
The LAPACK routine ZHEEV is used to diagonalise Hermitian Hamiltonians and obtain eigenvalues (band energies) and eigenvectors. The code then outputs the energy dispersion, orbital projections, and RGB-coloured band plots for visualisation.

The program is written in Fortran (F90) and requires LAPACK and BLAS libraries. Outputs are intended for visualisation with Gnuplot or Matplotlib.

Required input files:
BiTeI_hr_triv.dat (trivial phase Hamiltonian)
BiTeI_hr_top.dat (topological phase Hamiltonian)
BiTeI.nnkp (reciprocal lattice vectors, in the “begin recip_lattice” section)

The workflow is as follows:

1. Read input Hamiltonians for trivial and topological phases.
2. Define a k-path between high-symmetry points (L → A → H) in the Brillouin zone.
3. Construct the k-space Hamiltonian via Fourier transform: `H(k) = Σ_R H(R) e^(i k·R).`
4. Interpolate between the two phases using `H_alpha = α H_top + (1 - α) H_triv`, where `α = 0.7754` simulates the critical intermediate phase.
5. Diagonalise each `H_alpha(k)` with ZHEEV to obtain eigenvalues and eigenvectors.
6. Calculate orbital projections for Bi, Te, and I orbitals, and write the results to several output files.

Output files:
bandtop.dat: band energies and orbital weights
bandrgbtop.dat: band energies with RGB colour coding
bandrgb.plt: Gnuplot script for plotting
bandrgb.pdf: generated band structure plot

To compile the program, ensure LAPACK and BLAS are installed and run:
`gfortran Wannier_band_structure.f90 -llapack -lblas -o Wannier_band_structure`

Then execute:
`./Wannier_band_structure`

To generate the band structure plot, run:
`gnuplot bandrgb.plt`



**Notes:** The prefix variable in the code (character(len=80):: prefix="BiTeI") can be changed to use a different material. The interpolation factor alpha controls how “topological” the phase is. The default k-path follows L → A → H in the hexagonal Brillouin zone. For spin-polarised or 3D band plots, the output can be extended and visualised with Python or Matplotlib.

This project simulates surface state effects and band inversion across a topological transition in BiTeI, using tight-binding models and LAPACK-based eigenvalue solvers. It was developed as part of an MPhys research project at the University of Manchester, focusing on computational modelling of topological insulators and Rashba systems.
