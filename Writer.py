import numpy as np

import Parameters

def writeCoords(step, coords, geomFile, rawFile):
    
    # Write XYZ format
    dimensions = Parameters.Parameters.dimensions
    geomFile.write( str(len(coords)) + "\n")
    geomFile.write( "Step: " + str(step) + "\n" )
    
    for n in range(len(coords)):
        outArray = np.array( [ coords[n,j] for j in range(len(coords[0])) ] )
        geomFile.write( "\t".join(map(str,np.round(outArray,5))) + "\n" )
    
    # For radial distribution function, write out raw data
    for n in range(len(coords)):
        
        outArray = np.array( [ coords[n,j] for j in range(len(coords[0])) ] )
        rawFile.write( "\t".join(map(str,np.round(outArray,5))) + "\n" )

def writeThermo(step, V, potFile):

    # Write potential as a function of time
    potFile.write( "%5.5e\t%5.5e\n" % (step,V) )


