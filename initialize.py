import numpy as np
import matplotlib.pyplot as plt
import Parameters

def initParticles( NParticles=Parameters.Parameters.NParticles, latticeType=Parameters.Parameters.latticeType, 
        dimensions=Parameters.Parameters.dimensions, latticeLength=Parameters.Parameters.latticeLength ):
    """
    Initialize particles according to specified lattice arrangement and size.

    INPUT:
    NParticles:     int     -- Number of total particles in simulation box (Default: 27)
                            -- Must be a square in 2D or a cube in 3D
    latticeType:    'SC'    -- Simple cubic lattice (Default: 'SC')
    latticeConst:   float   -- Length of unit cell (Dafault: 1.0)

    OUTPUT:
    coords:         nd.array() -- shape = NParticles x 4 = Type, x,y,z

    """
    L = latticeLength
    radius = Parameters.Parameters.particleDiameter / 2

    coords = []
    if(latticeType=="SC"):
        if ( dimensions == 2 ):
            # Get important properties
            L12 = L / NParticles ** (1/2)
            for i in np.arange( 0, L, L12 ):
                for j in np.arange( 0, L, L12 ):
                    coords.append([ 1,i,j ])
            n = len(coords) / L ** 2 # Particle density
            phi = len(coords) * np.pi * radius**2 / L**2 # Volume of particles / Volume of box
        
        elif ( dimensions == 3 ):
         # Get important properties
            L13 = L / NParticles ** (1/3)
            for i in np.arange( 0, L, L13 ):
                for j in np.arange( 0, L, L13 ):
                    for k in np.arange( 0, L, L13 ):
                        coords.append([ 1,i,j,k ])
            for i in np.arange( 0, L-L13, L13 ):
                for j in np.arange( 0, L-L13, L13 ):
                    for k in np.arange( 0, L-L13, L13 ):
                        coords.append([ 1,i+L13/2,j+L13/2,k+L13/2 ])
            n = len(coords) / L ** 3 # Particle density
            phi = len(coords) * (4/3) * np.pi * radius**3 / L**3
        coords = np.array(coords)
    elif(latticeType=="BCC"):
        L13 = L / NParticles ** (1/3)
        for i in np.arange( 0, L, L13 ): #edge atoms
            for j in np.arange( 0, L, L13 ):
                for k in np.arange( 0, L, L13 ):
                    coords.append([ 1,i,j,k ])
        for i in np.arange( 0, L-L13, L13 ): #center atoms
            for j in np.arange( 0, L-L13, L13 ):
                for k in np.arange( 0, L-L13, L13 ):
                    coords.append([ 1,i+L13/2,j+L13/2,k+L13/2 ])
        coords = np.array(coords)
        n = len(coords) / L ** 3 # Particle density
        phi = len(coords) * (4/3) * np.pi * radius**3 / L**3
    elif(latticeType=="FCC"):
        # Get important properties
        L13 = L / NParticles ** (1/3)
        for i in np.arange( 0, L, L13 ): #lattice edges atoms
            for j in np.arange( 0, L, L13 ):
                for k in np.arange( 0, L, L13 ):
                    coords.append([ 1,i,j,k ])
        for i in np.arange( 0, L-L13/2, L13 ): #top bottom atoms
            for j in np.arange( 0, L-L13/2, L13 ):
                for k in np.arange( 0, L, L13 ):
                    #coords.append([ 1,i*L13+L13/2,j*L13+L13/2,k*L13+L13/2 ])
                    coords.append([ 1,i+L13/2,j+L13/2,k ])
        for i in np.arange( 0, L, L13 ): #front rear atoms
            for j in np.arange( 0, L-L13/2, L13 ):
                for k in np.arange( 0, L-L13/2, L13 ):
                    #coords.append([ 1,i*L13+L13/2,j*L13+L13/2,k*L13+L13/2 ])
                    coords.append([ 1,i,j+L13/2,k +L13/2])
        for i in np.arange( 0, L-L13/2, L13 ): #left right atoms
            for j in np.arange( 0, L, L13 ):
                for k in np.arange( 0, L-L13/2, L13 ):
                    #coords.append([ 1,i*L13+L13/2,j*L13+L13/2,k*L13+L13/2 ])
                    coords.append([ 1,i+L13/2,j,k +L13/2])
        for i in np.arange( 0, L-L13, L13/2 ): #center atoms
            for j in np.arange( 0, L-L13, L13/2 ):
                for k in np.arange( 0, L-L13, L13/2 ):
                    #coords.append([ 1,i*L13+L13/2,j*L13+L13/2,k*L13+L13/2 ])
                    coords.append([ 1,i+L13/2,j+L13/2,k+L13/2 ])
        n = len(coords) / L ** 3 # Particle density
        phi = len(coords) * (4/3) * np.pi * radius**3 / L**3
        coords = np.array(coords)
     

    #print ("The particle density is %5.4f" % (phi) )
    np.savetxt("geometry0.xyz", coords)

    print ("\n\tNumber Density:", n)
    print ("\n\tVolume Fraction:", phi)

    #print ("The length of coordinate array is:", len(coords))
    # Define boundary points
    Lmin = np.array([ np.min(coords[:,d+1]) for d in range(dimensions) ]) - Parameters.Parameters.particleDiameter
    Lmax = np.array([ np.max(coords[:,d+1]) for d in range(dimensions) ]) + Parameters.Parameters.particleDiameter
    print ("Max and Min")
    print (Lmax)
    print (Lmin)
    return coords, Lmin, Lmax 


if ( __name__ == "__main__" ):
    initParticles()
  
