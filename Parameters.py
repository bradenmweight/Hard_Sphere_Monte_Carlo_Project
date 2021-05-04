import numpy as np

class Parameters():

    dimensions = 3
    NParticles = 100
    particleDiameter = 1.0
    latticeType="SC" # Only "SC", "FCC" is implemented
    #latticeLength=9
    phi0 = 0.20
    NSteps = int(10**6)
    stepSize = 0.2
    kb, Temp = 1, 1 # Reduced Units

    boundaryType = 'PBC' # "Fixed", "PBC"
    potentialType = "LJ" # "Hard" , "LJ", "Free"

    # For LJ potential, define well parameter
    # Sigma taken to be diameter of particles
    eps = 0.5




