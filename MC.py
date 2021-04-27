import numpy as np
import random

import Parameters
from Writer import writeCoords
from Potential import getProb, getProbSingleParticle

def getGlobals():
    global dimensions, NSteps, stepSize, Diam, boundaryType
    dimensions = Parameters.Parameters.dimensions
    NSteps = Parameters.Parameters.NSteps
    stepSize = Parameters.Parameters.stepSize
    Diam = Parameters.Parameters.particleDiameter
    boundaryType = Parameters.Parameters.boundaryType

    # Initialize Output Files
    global geomFile, rawFile
    geomFile = open("dynamics.xyz","w")
    rawFile = open("dynamics.raw","w")


def getStep( coordsOLD, Lmin, Lmax, step, POLD ):

    # STEP 1
    # Get uniform random numbers for each particle in each dimension
    #   and move either forward/backward/left/right/up/down

    coordsNEW = coordsOLD * 1

    for n in range( len(coordsNEW) ): # Loop over particles
        for d in range( dimensions ):
            #if ( random.random() < 0.5 ):
            coordsNEW[n,d+1] += (random.random()*2-1) * stepSize# * random.random() # Here I am implenting non-uniform step size
            #else:
            #    coordsNEW[n,d+1] -= stepSize# * random.random()
            
            if ( boundaryType == "PBC" ):
                coordsNEW[n,d+1] = coordsNEW[n,d+1] % (Lmax[d]-Lmin[d])

    

    # STEP 2
    # Compute probability ratio of current and next step

    PNEW = getProb(coordsNEW, Lmin, Lmax)
    ratio = PNEW / POLD

    # STEP 3
    # Check to see if we accept new step based on ratio of probabilities

    if ( random.random() < ratio ):
        return coordsNEW, PNEW
    else:
        print ("SKIPPED.")
        return coordsOLD, POLD

def getStepSingleParticle( coordsOLD, Lmin, Lmax, step, POLD, ind ):

    # STEP 1
    # Get uniform random numbers for each particle in each dimension
    #   and move either forward/backward/left/right/up/down

    coordsNEW = coordsOLD * 1

    for d in range( dimensions ):
        coordsNEW[ind,d+1] += (random.random()*2-1) * stepSize
        
        if ( boundaryType == "PBC" ):
            coordsNEW[ind,d+1] = coordsNEW[ind,d+1] % (Lmax[d]-Lmin[d])

    

    # STEP 2
    # Compute probability ratio of current and next step

    PNEW = getProbSingleParticle( coordsNEW, Lmin, Lmax, ind )
    ratio = PNEW / POLD

    # STEP 3
    # Check to see if we accept new step based on ratio of probabilities

    if ( random.random() < ratio ):
        return coordsNEW, PNEW
    else:
        #print ("SKIPPED.")
        return coordsOLD, POLD

def runMCSingleParticleVersion( coords, Lmin, Lmax ):
    """
    This runs the Monte Carlo propagation.
    """

    getGlobals()

    print ("\nStarting MCMC algorithm with %5.0f steps and %5.0f particles" %(NSteps, len(coords)) )

    POLD = getProb(coords, Lmin, Lmax)

    for step in range(NSteps):
        print ("Step:", step)
        writeCoords(step,coords,geomFile,rawFile)
        for n in range(len(coords)):
            coords, POLD = getStepSingleParticle( coords, Lmin, Lmax, step, POLD, n )

    return None

def runMC( coords, Lmin, Lmax ):
    """
    This runs the Monte Carlo propagation.
    """

    getGlobals()

    print ("\nStarting MCMC algorithm with %5.0f steps and %5.0f particles" %(NSteps, len(coords)) )

    POLD = getProb(coords, Lmin, Lmax)

    for step in range(NSteps):
        print ("Step:", step)
        writeCoords(step,coords,geomFile,rawFile)
        coords, POLD = getStep( coords, Lmin, Lmax, step, POLD )

    return None


if ( __name__ == "__main__" ):
    coords = np.array([ [0,0,0], [0,0,1], [1,1,1] ])
    getMC( coords )
