# MATLAB_ED_BHM

This project compiles some functions written in MATLAB and are used for exact diagonalization simulation of 1D, bose-hubbard systems to compare with exeperiments of bosonic atoms in optical lattices in quantum gas microscopes. Note that this then implies all simulations performed by the basic functions in this code are inherently for closed systems since they will have only unitary dynamics.

## Contained Functions/Uses

[ ] Add n-point correlation functions to be called

[ ] Add monte-carlo sampling for fast and comparable data sampling

### BasisMake.m

This efficiently constructs the bose-hubbard basis states for a given number of particles "NPart" and number of lattice sites "NSites"

### cabFx.m

This calcultes the configurational correlation function for a  many-body state in a bose-hubbard system projected onto the Fock state basis. For details of physical use and importance look in SI of https://arxiv.org/abs/1805.09819

### DensityEval.m

Returns the on-site density average and on-site density variance from the wave function evaluated.

### ExactDiagTimeFx.m

Input of the initial state, evaluation times, and hamiltonian to compute exact diagonalization for efficient calculation of long-time dynamics and eigenstates. Output is time evolved wave function, eigenstates, and eigenenergies.

### FindTunnelingTerms.m

Deprecated function

### HilbDim.m

Returns hilbert space dimension from number of particles (NPart) and sites (NSites) in the bose-hubbard hamiltonian

### MakeDisorderHam.m

Constructs the disorder hamiltonian given different desired disorder distributions and basis states.

### MakeHamiltoniansAndBasis.m

If it doesn't already exist, this computes the interaction and tunneling hamiltonian matrices for a given bose-hubbard 1-D system given number of lattice sites, number of particles, boundary conditions, and tunneling (nearest neighbor, next-nearest neighbor, etc.).


## Example Script for single-particle physics

### andersonLocalizationScript.m

Uses exact diagonalization functions to compare single-particle localization in the bose-hubbard model for varius disorders strength and types.

## manybodyLocalizationScript.m

Uses exact diagonalization functions to compare many-body localization in the bose-hubbard model for varius disorders strengths, types, and interactions. 

[ ] Add larger system size comparisons

[ ] Add many-body to anderson localization comparisons


## Hamiltonians

This simply saves all previously constructed hamiltonians to save time for future exact diagonalization calculations
