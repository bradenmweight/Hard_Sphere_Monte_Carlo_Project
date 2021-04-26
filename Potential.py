import numpy as np

import Parameters

def getV( coords, Lmin, Lmax ):
    dimensions = Parameters.Parameters.dimensions
    V = 0

    """
    # Non-interacting
    kbT = Parameters.Parameters.kb * Parameters.Parameters.Temp
    V = (3/2) * Parameters.Parameters.NParticles * Parameters.Parameters.Temp
    """

    # Hard Sphere Potential
    Diam = Parameters.Parameters.particleDiameter
    for n in range( len(coords) ):
        for m in range( n+1, len(coords) ):
            dR = np.array([ coords[n,d+1] - coords[m,d+1] for d in range(dimensions) ])
            r = np.linalg.norm( dR )
            if ( r < Diam ):
                V += 10**8

    if ( Parameters.Parameters.boundaryType == "Fixed" ):
        # Fixed boundary walls
        for n in range( len(coords) ):
            # Bottom
            dR = np.array([ coords[n,d+1] - Lmin[d] for d in range(dimensions) ])
            r = np.linalg.norm( dR )
            if ( r < Parameters.Parameters.particleDiameter / 2 or dR.any() < 0 ):
                V += 10**8
                print ("WALL BOTTOM")
            # Top
            dR = np.array([ Lmax[d] - coords[n,d+1] for d in range(dimensions) ])
            r = np.linalg.norm( dR )
            if ( r < Parameters.Parameters.particleDiameter / 2 or dR.any() < 0 ):
                print ("WALL TOP")
                V += 10**8
            
    return V

def getProb( coords, Lmin, Lmax ):
    V = getV( coords, Lmin, Lmax )
    kb = Parameters.Parameters.kb
    T = Parameters.Parameters.Temp
    return np.exp( - V / kb / T)

def getLogP( coords, Lmin, Lmax ):
    V = getV( coords, Lmin, Lmax )
    kb = Parameters.Parameters.kb
    T = Parameters.Parameters.Temp
    return - V / kb / T


if ( __name__ == "__main__" ):

    coords = np.array([ [1.0,0,0,0], [1.0, 0, 0, 2], [1.0, 1,3,1] ])
    print ("Natoms =", len(coords))
    V = getV( coords )
    P = getProb( coords )
    LogP = getLogP( coords )
    print ("V =", V)
    print ("P =", P)
    print ("LogP =", LogP)
