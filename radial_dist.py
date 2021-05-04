import numpy as np
import sys

import Parameters


def getRDFHist( coords, L, NBINS, BINS, boundaryType ):
    """
    Return the normalized RDF for a single snapshot.
    """

    dim = Parameters.Parameters.dimensions
    diameter = Parameters.Parameters.particleDiameter

    gr = np.zeros(( NBINS ))
    for i in range( len(coords) ):
        for j in range( i+1, len(coords) ):
            if ( boundaryType == "PBC" ):
                dR = np.zeros(dim)
                for d in range(dim):
                    diff = coords[i,d+1] - coords[j,d+1]
                    if diff > L/2:
                        dR[d] = L - diff
                    elif diff <= -L/2:
                        dR[d] = L + diff
                    else: 
                        dR[d] = diff
            else:
                dR = np.array([ coords[i,d+1] - coords[j,d+1] for d in range(dim) ])
            r = np.linalg.norm( dR )
            for b in range( NBINS-1 ):
                if ( BINS[b] < r and r < BINS[b+1] ):
                    gr[b] += 1
    
    # Normalize histogram by spherical shell*dr and non-interacting density

    for b in range(NBINS-1):
        #Adr = 4 * np.pi * BINS[b] ** 2 * dr * ( dim == 3 ) +  2 * np.pi * BINS[b] * dr * ( dim == 2 )
        Adr = (4/3) * np.pi * (BINS[b+1] ** 3 - BINS[b] ** 3) * ( dim == 3 ) +  np.pi * (BINS[b+1] ** 2 - BINS[b] ** 2) * ( dim == 2 )
        gr[b] /= Adr

    return gr

if ( __name__ == "__main__" ):
    dimensions = Parameters.Parameters.dimensions
    NSteps = Parameters.Parameters.NSteps
    diameter = Parameters.Parameters.particleDiameter
    boundaryType = Parameters.Parameters.boundaryType
    #L = Parameters.Parameters.latticeLength

    NSteps = 50000//100

    if ( len(sys.argv) == 2 and len(sys.argv[1].split(".")) < 2 ):
        name = sys.argv[1]
        print (f"Reading File: {name}.xyz AND {name}.raw")
    else:
        name = "dynamics"

    NAtoms = int( open(f"{name}.xyz","r").readlines()[0] )
    coords = np.loadtxt(f"{name}.raw").flatten().reshape((NSteps,NAtoms,dimensions+1))
   

    L = 8.05996  # THIS NEEDS TO BE CHANGED EACH TIME !!!!!!!!!! BE CAREFUL !!!!!!!!!!
    print ("\n\tBox Size (L):",L)

    dBIN = 0.1
    BINS = np.arange( 0, L/2, dBIN ) # I have assumed cubic box
    NBINS = len(BINS)
    gr = np.zeros(( len(BINS) ))
    
    NStart = 1000
    NSkip = 500
    for step in range( NStart,NSteps ):
        if ( step % NSkip == 0 ):
            print (f"Step: {step} of {NSteps}")
            gr += getRDFHist( coords[step], L, NBINS, BINS, boundaryType )
    
    NFrames = (NSteps - NStart) / NSkip
    gr /= NFrames

    V = L ** 2 * ( dimensions == 2) + L ** 3 * ( dimensions == 3 )
    gr *= V/NAtoms**2

    gr *= 2 # This is an ad-hoc correction

    np.savetxt("RDF.dat", np.array([BINS[:-1], gr[:-1]]).T )


