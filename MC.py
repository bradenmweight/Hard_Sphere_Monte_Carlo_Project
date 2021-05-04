import numpy as np
import random

import Parameters
from Writer import writeCoords, writeThermo
from Potential import getProb, getV

def getGlobals():
    global dimensions, NSteps, stepSize, Diam, boundaryType
    dimensions = Parameters.Parameters.dimensions
    NSteps = Parameters.Parameters.NSteps
    stepSize = Parameters.Parameters.stepSize
    Diam = Parameters.Parameters.particleDiameter
    boundaryType = Parameters.Parameters.boundaryType
    #L = Parameters.Parameters.latticeLength

    # Initialize Output Files
    global geomFile, rawFile, potFile
    geomFile = open("dynamics.xyz","w")
    rawFile = open("dynamics.raw","w")
    potFile = open("thermo.dat","w")


def getDists(coords):
    distx = np.subtract.outer(coords[:,1],coords[:,1])
    disty = np.subtract.outer(coords[:,2],coords[:,2])
    distz = np.subtract.outer(coords[:,3],coords[:,3])

    distr = np.sqrt( distx**2 + disty**2 + distz**2 )

    return distr


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

def getStepSingleParticle( coordsOLD, Lmin, Lmax, step, VOLD ):

    coordsNEW = coordsOLD * 1
    L = Lmax - Lmin

    a = random.randint(0,len(coordsNEW)-1) # Choose random particle to move
    for d in range( dimensions ):
        coordsNEW[a,d+1] += (random.random()*2-1) * stepSize
        
        if ( boundaryType == "PBC" ):
            coordsNEW[a,d+1] = coordsNEW[a,d+1] % L[d]

    VNEW = getV( coordsNEW, Lmin, Lmax )

    if ( VNEW < VOLD ):
        return coordsNEW, VNEW
    else:
        kb = Parameters.Parameters.kb
        T = Parameters.Parameters.Temp
        P = np.exp(-(VNEW-VOLD)/kb/T)
        if ( random.random() < P ):
            return coordsNEW, VNEW
        else:
            print ("SKIPPED.")
            return coordsOLD, VOLD

       

def runMCSingleParticleVersion( coords, Lmin, Lmax ):
    """
    This runs the Monte Carlo propagation.
    """

    getGlobals()
    T = Parameters.Parameters.Temp
    kb = Parameters.Parameters.kb

    print ("\nStarting MCMC algorithm with %5.0f steps and %5.0f particles" %(NSteps, len(coords)) )

    VOLD = 1.0

    for step in range(NSteps):
        print ("Step:", step)
        if ( step % 100 == 0):
            writeCoords(step,coords,geomFile,rawFile)
        if ( step > 0 ):
            #print ("V = %5.4f" % (-np.log(POLD)*kb*T) )
            writeThermo( step, VOLD, potFile ) # Write poential energy at each full step -- after all particles have chance to move
        coords, VOLD = getStepSingleParticle( coords, Lmin, Lmax, step, VOLD )

    return None

def runMC( coords, Lmin, Lmax ):
    """
    This runs the Monte Carlo propagation.
    """

    getGlobals()
    T = Parameters.Parameters.Temp
    kb = Parameters.Parameters.kb

    print ("\nStarting MCMC algorithm with %5.0f steps and %5.0f particles" %(NSteps, len(coords)) )

    POLD = getProb(coords, Lmin, Lmax)

    for step in range(NSteps):
        print ("Step:", step)
        writeCoords(step,coords,geomFile,rawFile)
        writeThermo( step, -np.log(POLD)*kb*T, potFile ) # Write poential energy at each full step -- after all particles have chance to move
    
        coords, POLD = getStep( coords, Lmin, Lmax, step, POLD )

    return None


if ( __name__ == "__main__" ):
    coords = np.array([ [0,0,0], [0,0,1], [1,1,1] ])
    getMC( coords )
