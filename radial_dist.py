import numpy as np

import Parameters


def getRDFHist( coords, Lmin, Lmax, NBINS, boundaryType ):
    """
    Given a dynamics trajectory, this will return the normalized RDF
    """
    diameter = Parameters.Parameters.particleDiameter
    boundaryType = Parameters.Parameters.boundaryType

    BINS = np.linspace( 0,(Lmax[0]-Lmin[0])/2,NBINS ) # I have assumed cubic box
    HIST = np.zeros(( NBINS ))
    for i in range( len(coords) ):
        for j in range( i+1, len(coords) ):
            if ( boundaryType == "PBC" ):
                dR = np.array([ coords[i,d+1] % (Lmax[d]-Lmin[d]) - coords[j,d+1] % (Lmax[d]-Lmin[d]) for d in range(dimensions) ])
            else:
                dR = np.array([ coords[i,d+1] - coords[j,d+1] for d in range(dimensions) ])
            r = np.linalg.norm( dR )
            for b in range( NBINS-1 ):
                if ( BINS[b] < r and r < BINS[b+1] ):
                    HIST[b] += 1
    
    # Normalize histogram by spherical shell*dr
    dr = BINS[1] - BINS[0]
    dim = Parameters.Parameters.dimensions
    L = Lmax[0]-Lmin[0]
    n = len(coords) / L ** ( 2 * ( dim == 2 ) + 3 * ( dim == 3 ) )
    for b in range(1,NBINS):
        Adr = 4 * np.pi * BINS[b] ** 2 * dr * ( dim == 3 ) +  2 * np.pi * BINS[b] * dr * ( dim == 2 )
        HIST[b] = HIST[b] / Adr

    return np.array([ BINS/diameter, HIST ])

if ( __name__ == "__main__" ):
    dimensions = Parameters.Parameters.dimensions
    NSteps = Parameters.Parameters.NSteps
    diameter = Parameters.Parameters.particleDiameter
    boundaryType = Parameters.Parameters.boundaryType

    NAtoms = int( open("dynamics.xyz","r").readlines()[0] )
    coords = np.loadtxt("dynamics.raw").flatten().reshape((NSteps,NAtoms,dimensions+1))

    Lmin = np.array([ np.min(coords[0,:,d+1]) for d in range(dimensions) ]) - Parameters.Parameters.particleDiameter
    Lmax = np.array([ np.max(coords[0,:,d+1]) for d in range(dimensions) ]) + Parameters.Parameters.particleDiameter
    
    print ("\n\tBox Size (L):",Lmax[0]-Lmin[0])

    NBINS = 100
    BINS = np.linspace( 0,(Lmax[0]-Lmin[0])/2,NBINS ) # I have assumed cubic box
    HIST = np.zeros(( 2,len(BINS) ))
    
    NStart = int(NSteps / 10)
    NSkip = 50
    for step in range( NStart,NSteps ):
        if ( step % NSkip == 0 ):
            print (f"Step: {step} of {NSteps}")
            coords_curr = coords[step]
            HIST += getRDFHist( coords_curr, Lmin, Lmax, NBINS, boundaryType )
    
    HIST = HIST / ( (NSteps - NStart) / NSkip )
    np.savetxt("RDF.dat",HIST.T )

    Npts = 4096
    SSF = np.abs( np.fft.fft( np.append(HIST[1,:-1],np.zeros(( Npts - len(HIST[1,:-1]) ))),norm="ortho") )
    q = np.fft.fftfreq( Npts ) * 2 * np.pi / (BINS[1]-BINS[0])
    np.savetxt("SSF.dat", np.array([q,SSF]).T )
