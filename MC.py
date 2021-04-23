import numpy as np
import random

import Parameters
from Writer import writeCoords
from Potential import getProb

def getGlobals():
    global dimensions, NSteps, stepSize
    dimensions = Parameters.Parameters.dimensions
    NSteps = Parameters.Parameters.NSteps
    stepSize = Parameters.Parameters.stepSize

    # Initialize Output Files
    global geomFile, rawFile
    geomFile = open("dynamics.xyz","w")
    rawFile = open("dynamics.raw","w")


def getStep( coordsOLD, Lmin, Lmax, step, POLD ):

    # If first step, get initial probability    
    if ( step == 0 ):
        POLD = getProb(coordsOLD, Lmin, Lmax)

    # STEP 1
    # Get uniform random numbers for each particle in each dimension
    #   and move either forward/backward/left/right/up/down

    coordsNEW = coordsOLD * 1

    for n in range( len(coordsNEW) ): # Loop over particles
        for d in range( dimensions ):
            if ( random.random() < 0.5 ):
                coordsNEW[n,d+1] += stepSize * random.random() # Here I am implenting non-uniform step size
            else:
                coordsNEW[n,d+1] -= stepSize * random.random()

    # STEP 2
    # Compute probability ratio of current and next step

    PNEW = getProb(coordsNEW, Lmin, Lmax)
    ratio = PNEW / POLD

    """
    # Use for debugging if encountering nan or inf.
    if ( np.isnan(ratio) ):
        print ("\t'ratio' is nan.  Killing.\n")
        exit()
    elif( np.isinf(ratio) ):
        print ("\t'ratio' is inf.  Killing.\n")
        exit()
    """

    # STEP 3
    # Check to see if we accept new step based on ratio of probabilities

    if ( random.random() < ratio ):
        return coordsNEW, PNEW
    else:
        print ("SKIPPED.")
        return coordsOLD, POLD



def runMC( coords, Lmin, Lmax ):
    """
    This runs the Monte Carlo propagation.
    """

    getGlobals()

    print ("\nStarting MCMC algorithm with %5.0f steps." %(NSteps) )

    Prob = 0

    for step in range(NSteps):
        print ("Step:", step)
        writeCoords(step,coords,geomFile,rawFile)
        coords, Prob = getStep( coords, Lmin, Lmax, step, Prob )  

    return None


if ( __name__ == "__main__" ):
    coords = np.array([ [0,0,0], [0,0,1], [1,1,1] ])
    getMC( coords )
