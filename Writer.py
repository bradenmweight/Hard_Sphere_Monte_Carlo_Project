import numpy as np

def writeCoords(step, coords, geomFile):
    geomFile.write( str(len(coords)) + "\n")
    geomFile.write( "Step: " + str(step) + "\n" )
    for n in range(len(coords)):
        outArray = np.array( [1.0, coords[n,0], coords[n,1], coords[n,2] ] )
        geomFile.write( "\t".join(map(str,np.round(outArray,5))) + "\n" )
    #geomFile.write( "\n" )








