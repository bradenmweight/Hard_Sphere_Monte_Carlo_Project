import numpy as np
import random as rand
import math

# Define Constants
T = 1 # K
k = 1 # 1.380649 x 10-23 JK-1.
N = 10
diameter = 0.1
# walls, particle interactions

def getProb(pos):
    E = np.ones(N) #energy doesnt change over time so one energy per particle
    for i in range(N-1):
        if ( abs(np.sqrt((pos[i,0])**2+(pos[i,1])**2+(pos[i,2])**2)-np.sqrt((pos[i+1,0])**2+(pos[i+1,1])**2+(pos[i+1,2])**2)) < diameter ):
            print("overlap!")
            V = math.inf
            E[i] += V
            #P = 0 
        else:
            V = 0 
            
    H = np.sum(E*N) #3/2KT would be constant since we're not changing T 
    P = np.exp(-E/k/T)/np.sum(np.exp(-E/k/T)) 
    
    return P 

#everything below here is for testing purposes:
"""
pos = []
for n in range(N):
    pos.append( [ rand.random(), rand.random(), rand.random() ] )    
pos = np.array(pos)

test = getProb(pos)
print(test)
"""
