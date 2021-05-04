import numpy as np
import emcee
import sys


def getData():

    data = np.loadtxt("thermo.dat")
    return data
    
def getAutoCorrelation(timeseries):

    print ("\n\tGetting Correlation Function")

    time = timeseries[:,0]
    V = timeseries[:,1]

    AUTO_V = emcee.autocorr.function_1d( V )
    tau_V = emcee.autocorr.integrated_time(V, c=5, tol=5, quiet=True)

    print ("\tCorrelation Time 'tau_V':", int(tau_V[0]), "of", int(len(time)) )

    return tau_V[0], np.array([ time, AUTO_V ]).T

def trim_Data(timeseries, tau_V):

    print ("\tI am slicing the data.")

    trimmed_data = timeseries[int(tau_V)::int(tau_V/2)]

    return trimmed_data


def printResults(trimmed_data, AUTO_V):

    print ("\tI am printing the results to file.")

    np.savetxt( "thermo_sliced.dat", trimmed_data )
    #np.savetxt("V_auto-correlation.dat", AUTO_V) # This is a large file
    np.savetxt( "V_auto-correlation_shortened.dat", AUTO_V[::100] )


if ( __name__ == "__main__" ):

    timeseries = getData()
    tau_V, AUTO_V = getAutoCorrelation(timeseries)
    trimmed_data = trim_Data(timeseries, tau_V)
    printResults(trimmed_data, AUTO_V)

    print ("\tHave a nice day! :)\n")
