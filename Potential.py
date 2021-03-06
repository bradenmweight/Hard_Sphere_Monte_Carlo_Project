import numpy as np

import Parameters

def getBoundaryV( coords, dimensions):
    # Fixed boundary walls
    for n in range( len(coords) ):
        # Bottom
        dR = np.array([ coords[n,d+1] - Lmin[d] for d in range(dimensions) ])
        r = np.linalg.norm( dR )
        if ( r < Parameters.Parameters.particleDiameter / 2 or dR.any() < 0 ):
            print ("WALL BOTTOM")
            V += 10**3
            return V
        # Top
        dR = np.array([ Lmax[d] - coords[n,d+1] for d in range(dimensions) ])
        r = np.linalg.norm( dR )
        if ( r < Parameters.Parameters.particleDiameter / 2 or dR.any() < 0 ):
            print ("WALL TOP")
            V += 10**3
            return V

def getV( coords, Lmin, Lmax ):
    dimensions = Parameters.Parameters.dimensions
    potentialType = Parameters.Parameters.potentialType
    boundaryType = Parameters.Parameters.boundaryType
    L = Lmax-Lmin
    V = 0

    ##### DEFINE WALL BOUNDARIES #####
    if ( Parameters.Parameters.boundaryType == "Fixed" ):
        V = getBoundaryV( coords, dimensions)

    if ( potentialType == "Free" ):
        # Non-interacting
        kbT = Parameters.Parameters.kb * Parameters.Parameters.Temp
        V = (dimensions/2) * Parameters.Parameters.NParticles * Parameters.Parameters.Temp
        return V

    elif ( potentialType == "Hard" ):
        # Hard Sphere/Disc Potential
        Diam = Parameters.Parameters.particleDiameter
        for n in range( len(coords) ):
            for m in range( n+1, len(coords) ):
                if ( boundaryType == "PBC" ):
                    #dR = np.array([ coords[m,d+1] % L[d] - coords[n,d+1] % (Lmax[d]-Lmin[d]) for d in range(dimensions) ])
                    dR = np.zeros(dimensions)
                    for d in range(dimensions):
                        diff = (coords[m,d+1] - coords[n,d+1])
                        if diff > L[d]/2:
                            dR[d] = L[d] - diff
                        elif diff <= -L[d]/2:
                            dR[d] = L[d] + diff
                        else: 
                            dR[d] = diff
                else:
                    dR = np.array([ coords[n,d+1] - coords[m,d+1] for d in range(dimensions) ])
                r = np.linalg.norm( dR )
                if ( r < Diam ):
                    V += 10**6
        return V

    elif( potentialType == "LJ" ):
        # Lennard-Jones Potential
        Diam = Parameters.Parameters.particleDiameter
        eps = Parameters.Parameters.eps
        for n in range( len(coords) ):
            for m in range( n+1, len(coords) ):
                if ( boundaryType == "PBC" ):
                    dR = np.zeros(dimensions)
                    for d in range(dimensions):
                        diff = (coords[m,d+1] - coords[n,d+1])
                        if diff > L[d]/2:
                            dR[d] = L[d] - diff
                        elif diff <= -L[d]/2:
                            dR[d] = L[d] + diff
                        else: 
                            dR[d] = diff
                else:
                    dR = np.array([ coords[n,d+1] - coords[m,d+1] for d in range(dimensions) ])
                r = np.linalg.norm( dR )
                V += 4*eps*( (Diam/r)**12 - (Diam/r)**6 )
        return V

    print ("I should not get here!")
    print ("Potential Type:", potentialType)
    exit()

def getVSingleParticle( coords, Lmin, Lmax, ind ):
    dimensions = Parameters.Parameters.dimensions
    potentialType = Parameters.Parameters.potentialType
    boundaryType = Parameters.Parameters.boundaryType
    L = Lmax-Lmin
    V = 0

    ##### DEFINE WALL BOUNDARIES #####
    if ( Parameters.Parameters.boundaryType == "Fixed" ):
        V = getBoundaryV( coords, dimensions)

    if ( potentialType == "Free" ):
        # Non-interacting
        kbT = Parameters.Parameters.kb * Parameters.Parameters.Temp
        V = (dimensions/2) * Parameters.Parameters.NParticles * Parameters.Parameters.Temp
        return V

    elif ( potentialType == "Hard" ):
        # Hard Sphere/Disc Potential
        Diam = Parameters.Parameters.particleDiameter
        for n in range( len(coords) ):
            if ( n != ind ):
                if ( boundaryType == "PBC" ):
                    dR = np.zeros(dimensions)
                    for d in range(dimensions):
                        diff = (coords[n,d+1] - coords[ind,d+1])
                        if diff > (Lmax[d]-Lmin[d])/2:
                            dR[d] = (Lmax[d]-Lmin[d]) - diff
                        elif diff <= -(Lmax[d]-Lmin[d])/2:
                            dR[d] = (Lmax[d]-Lmin[d]) + diff
                        else: 
                            dR[d] = diff
                else:
                    dR = np.array([ coords[n,d+1] - coords[m,d+1] for d in range(dimensions) ])
                r = np.linalg.norm( dR )
                if ( r < Diam ):
                    V += 10**10
        return V

    elif( potentialType == "LJ" ):
        # Lennard-Jones Potential
        Diam = Parameters.Parameters.particleDiameter
        sig = Diam
        eps = Parameters.Parameters.eps
        for n in range( len(coords) ):
            if ( n != ind ):
                if ( boundaryType == "PBC" ):
                    dR = np.zeros(dimensions)
                
                    for d in range(dimensions):
                    #    dR[d] = ( coords[n,d+1] - coords[ind,d+1] ) % (L[d]/2)
                        #print (dR[d])
                        diff = coords[n,d+1] - coords[ind,d+1]
                        if diff > (Lmax[d]-Lmin[d])/2:
                            dR[d] = (Lmax[d]-Lmin[d]) - diff
                        elif diff <= -(Lmax[d]-Lmin[d])/2:
                            dR[d] = (Lmax[d]-Lmin[d]) + diff
                        else: 
                            dR[d] = diff
               
                r = np.linalg.norm( dR )
                V += 4*eps*( (sig/r)**12 - (sig/r)**6 )
        return V

def getProb( coords, Lmin, Lmax ):
    V = getV( coords, Lmin, Lmax )
    kb = Parameters.Parameters.kb
    T = Parameters.Parameters.Temp
    return np.exp( - V / kb / T)

if ( __name__ == "__main__" ):

    coords = np.array([ [1.0,0,0,0], [1.0, 0, 0, 2], [1.0, 1,3,1] ])
    print ("Natoms =", len(coords))
    V = getV( coords )
    P = getProb( coords )
    LogP = getLogP( coords )
    print ("V =", V)
    print ("P =", P)
    print ("LogP =", LogP)
