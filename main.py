# This will be the main file for this package.
# Chloe can edit

import numpy as np
import initialize




# Initialize Particles/Box Somehow: Lattice or random? Periodic boundaries or walls?
# Write Hamiltonian and compute initial energy of system.

# For step in NSteps:
    # Take a Monte Carlo step.
    # Calculate new energy.
    # Compute ratio of Boltzmann probabilities
    # Accept step if ratio > random [0,1]
    # If accepted, compute necessary information: pressure, temperature, etc. g(r), SSF, etc.

class Parameters():

    dimensions = 2
    NParticles = 22
    particleDiameter = 1.0
    latticeType="SC"
    latticeLength=20



if ( __name__ == "__main__" ):

    coords = initialize.initParticles( Parameters.NParticles, Parameters.latticeType, Parameters.dimensions, Parameters.latticeLength )












