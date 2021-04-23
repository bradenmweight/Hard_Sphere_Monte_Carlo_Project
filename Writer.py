import numpy as np

def writeCoords(step, coords, geomFile):
    for n in range(len(coords)):
        outArray = [step, coords[n,0], coords[n,1], coords[n,2] ]
        geomFile.write( "\t".join(map(str,outArray)) + "\n" )


