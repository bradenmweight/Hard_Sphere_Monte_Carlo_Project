import numpy as np

class Parameters():

    ### SIMULATION DETAILS ###
    dimensions = 2
    NParticles = 100 # Code will round up to nearest square (2D) or cube (3D)
    particleDiameter = 1.0
    latticeType="SC" # Only "SC" is implemented
    latticeLength=20
    NSteps = 5000
    stepSize = 0.2
    kb, Temp = 1, 2 # Reduced Units

    boundaryType = 'PBC' # "Fixed", "PBC"
    potentialType = "LJ" # "Hard" , "LJ", "Free"

    # For LJ potential, define well parameter
    # Sigma taken to be diameter of particles
    eps = 5.0



