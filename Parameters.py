import numpy as np

class Parameters():

    ### SIMULATION DETAILS ###
    dimensions = 3
    NParticles = 20 # Code will round up to nearest square or cube
    particleDiameter = 2.0
    latticeType="SC" # Only "SC" is implemented
    latticeLength=10
    NSteps = 10000
    stepSize = 0.2
    kb = 1
    Temp = 2 # Reduced Units

    boundaryType = 'Fixed' # "PBC" not yet implemented


