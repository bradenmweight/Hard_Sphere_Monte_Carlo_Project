# This will be the main file for this package.
# Chloe can edit

import numpy as np




# Initialize Particles/Box Somehow: Lattice or random? Periodic boundaries or walls?
# Write Hamiltonian and compute initial energy of system.

# For step in NSteps:
    # Take a Monte Carlo step.
    # Calculate new energy.
    # Compute ratio of Boltzmann probabilities
    # Accept step if ratio > random [0,1]
    # If accepted, compute necessary information: pressure, temperature, etc. g(r), SSF, etc.

def initParticles( NParticles=10, latticeType="SC", latticeConst=1.0 ):

    # For now, choose to initialize particles in simple cubic lattice. Add other options later.
    
    # Simple Cubic Lattice
    a = latticeConst
    coords = np.array()

    return None




if ( __name__ == "__main__" ):

    initParticles()







