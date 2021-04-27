# This will be the main file for this package.
# Chloe can edit

import numpy as np

import initialize
import MC
import Parameters



# Initialize Particles/Box Somehow: Lattice or random? Periodic boundaries or walls?
# Write Hamiltonian and compute initial energy of system.

# For step in NSteps:
    # Take a Monte Carlo step.
    # Calculate new energy.
    # Compute ratio of Boltzmann probabilities
    # Accept step if ratio > random [0,1]
    # If accepted, compute necessary information: pressure, temperature, etc. g(r), SSF, etc.

# ANALYSIS
#   1. Radial Distrution Function   --> radial_dist.dat
#   2. Static Structure Factor      --> 
#   3. Compute time-correlation    --> 
#   4. SSF --> Phase Diagram via Hansen and Verlet freezing criterion (height of first peak in SSF)

# POSSIBLE FUTURE IMPROVEMENTS
#   1. Implement periodic boundary conditions: boundaryType = "Fixed", "PBC"
#   2. Implement different initial lattices: lattiveType = "SC", "FCC", "BCC", "Diamond"
#   3. Additional potential types: Lennard-Jones (6-12)
#   4. Implement 1-Dimensional code

if ( __name__ == "__main__" ):

    # Initialize the particles
    coords, Lmin, Lmax  = initialize.initParticles()
    #MC.runMC( coords, Lmin, Lmax ) # This runs as follows: (I) Move all particles randomly (II) Compute probability
    MC.runMCSingleParticleVersion( coords, Lmin, Lmax )  # This runs as follows: (I) Move a single particle randomly (II) Compute probability




