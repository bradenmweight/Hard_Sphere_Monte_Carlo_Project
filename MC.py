import numpy as np
import random
import Parameters

from Writer import writeCoords

def getGlobals():
    global dimensions, NSteps, stepSize
    dimensions = Parameters.Parameters.dimensions
    NSteps = Parameters.Parameters.NSteps
    stepSize = Parameters.Parameters.stepSize

    # Initialize Output Files
    global geomFile
    geomFile = open("dynamics.xyz","w")


def getStep( coordsOLD, step, POLD ):

    # If first step, get initial probability    
    if ( step == 0 ):
        POLD = 1 #getProb(coords)


    # STEP 1
    # Get uniform random numbers for each particle in each dimension
    #   and move either forward/backward/left/right/up/down

    coordsNEW = coordsOLD * 1

    for n in range( len(coordsNEW) ): # Loop over particles
        for d in range( dimensions ):
            if ( random.random() < 0.5 ):
                coordsNEW[n,d] += stepSize
            else:
                coordsNEW[n,d] -= stepSize

    # STEP 2
    # Compute probability ratio of current and next step

    PNEW = 1 #getProb(coords)
    ratio = PNEW / POLD

    if ( np.isnan(ratio) ):
        print ("\t'ratio' is nan.  Killing.\n")
        exit()
    elif( np.isinf(ratio) ):
        print ("\t'ratio' is inf.  Killing.\n")
        exit()

    # STEP 3
    # Check to see if we accept new step based on ratio of probabilities

    if ( random.random() < ratio ):
        return coordsNEW, PNEW
    else:
        return coordsOLD, POLD



def runMC( coords ):
    """
    This runs the Monte Carlo propagation.
    """

    getGlobals()

    print ("\nStarting MCMC algorithm with %5.0f steps." %(NSteps) )

    Prob = 0

    for step in range(NSteps):
        print ("Step:", step)
        coords, Prob = getStep( coords, step, Prob )
        writeCoords(step,coords,geomFile)

    

    return None


if ( __name__ == "__main__" ):
    coords = np.array([ [0,0,0], [0,0,1], [1,1,1] ])
    getMC( coords )
