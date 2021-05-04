import numpy as np
from scipy.special import jn
import sys
from time import time

import Parameters

def getSSF_SLOW( coords, L, Nq, q, boundaryType ):
        dimensions = Parameters.Parameters.dimensions
        #ssf = np.zeros(( len(q) ), dtype = complex)
        ssf = np.zeros(( len(q) ))
        for qi in range (len(q)):
            for i in range(len(coords)):
                    print ("Atom i,",i, "of", len(coords))
                    for j in range(i+1, len(coords)):
                            dR = np.zeros((dimensions))
                            for d in range(dimensions):
                                diff = coords[i,d+1] - coords[j,d+1]
                                if diff > L/2:
                                    dR[d] = L - diff
                                elif diff <= -L/2:
                                    dR[d] = L + diff
                                else:
                                    dR[d] = diff
                            r = np.linalg.norm( dR )
                            #ssf[qi] += np.exp(1j*q[qi]*r)
                            ssf[qi] += np.sin( q[qi]*r ) / (q[qi] * r) * ( dimensions == 3 ) 
                            ssf[qi] += jn( 0, q[qi]*r ) / (q[qi] * r) * ( dimensions == 2 ) # This might not be right. Unsure of (q*r)^-1 prefactor.
                            
        return ssf

def getSSF( coords, L, Nq, q, boundaryType ): # Faster implementation by putting loop over q inside loop over particles
        dimensions = Parameters.Parameters.dimensions
        #ssf = np.zeros(( len(q) ), dtype = complex)
        ssf = np.zeros(( len(q) ))
        
        for i in range(len(coords)):
            print ("Atom i,",i, "of", len(coords))
            for j in range(i+1, len(coords)):
                    dR = np.zeros((dimensions))
                    for d in range(dimensions):
                        diff = coords[i,d+1] - coords[j,d+1]
                        if diff > L/2:
                            dR[d] = L - diff
                        elif diff <= -L/2:
                            dR[d] = L + diff
                        else:
                            dR[d] = diff
                    r = np.linalg.norm( dR )
                    ssf += np.sin( q*r ) / (q * r) * ( dimensions == 3 ) 
                    #for qi in range (len(q)):
                    #    #ssf[qi] += np.exp(1j*q[qi]*r)
                    #    ssf[qi] += np.sin( q[qi]*r ) / (q[qi] * r) * ( dimensions == 3 ) 
                    #    ssf[qi] += jn( 0, q[qi]*r ) / (q[qi] * r) * ( dimensions == 2 ) # This might not be right. Unsure of (q*r)^-1 prefactor.
                            
        return ssf

if ( __name__ == "__main__" ):
    dimensions = Parameters.Parameters.dimensions
    NSteps = Parameters.Parameters.NSteps
    diameter = Parameters.Parameters.particleDiameter
    boundaryType = Parameters.Parameters.boundaryType
    #L = Parameters.Parameters.latticeLength

    if ( len(sys.argv) == 2 and len(sys.argv[1].split(".")) < 2 ):
        name = sys.argv[1]
    else:
        name = "dynamics"
    print (f"Reading File: {name}.xyz AND {name}.raw")


    NSteps = 50000//100

    NAtoms = int( open(f"{name}.xyz","r").readlines()[0] )
    coords = np.loadtxt(f"{name}.raw").flatten().reshape((NSteps,NAtoms,dimensions+1))

    Nq = 200
    dq = 0.1
    q = np.arange( 0,5,(5)/Nq ) + dq # I have assumed cubic box
    #ssf = np.zeros(( len(q) ),dtype = complex)
    ssf = np.zeros(( len(q) ))

    L = 11.02587 # THIS NEEDS TO BE CHANGED EACH TIME !!!!!!!!!! BE CAREFUL !!!!!!!!!!
    print ("\n\tBox Size (L):",L)
    print (f"\t(qmin,qmax) = ({q[0]},{q[-1]})")
    print (f"\tdq = {q[1]-q[0]}")

    NStart = int( NSteps / 10)
    NSkip = 25
    for step in range( NStart,NSteps ):
        if ( step % NSkip == 0 ):
            print (f"Step: {step} of {NSteps}")
            start = time()
            ssf += getSSF( coords[step], L, Nq, q, boundaryType )
            print ("Time:", time() - start)

    ssf /=  (NSteps - NStart) / NSkip
    ssf /= NAtoms

    ssf += 1

    np.savetxt("SSF_sin_FAST.dat",np.array([q[1:],ssf[1:]]).T )

    """
    ssfRE = np.real(ssf)
    ssfIM = np.imag(ssf)
    ssfABS = np.sqrt(np.real(np.conj(ssf)*ssf))
    np.savetxt("SSF.dat",np.array([q,ssfRE,ssfIM,ssfABS]).T )
    """
