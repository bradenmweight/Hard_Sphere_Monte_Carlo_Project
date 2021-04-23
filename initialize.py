import numpy as np

def initParticles( NParticles=28, latticeType="SC", dimensions=3, latticeLength=10 ):
    """
    Initialize particles according to specified lattice arrangement and size.

    INPUT:
    NParticles:     int     -- Number of total particles in simulation box (Default: 27)
                            -- Must be a square in 2D or a cube in 3D
    latticeType:    'SC'    -- Simple cubic lattice (Default: 'SC')
    latticeConst:   float   -- Length of unit cell (Dafault: 1.0)

    OUTPUT:
    coords:         nd.array()

    """
    
    # Simple Cubic Lattice
    assert ( latticeType == "SC" and dimensions != 1 ), "Only 2D and 3D SC lattice is currently implemented."
    
    L = latticeLength

    coords = []
    if ( dimensions == 2 ):
        # Get important properties
        L12 = L / NParticles ** (1/2)
        phi = NParticles / L ** 2 # Particle density

        for i in np.arange( 0, L, L12 ):
            for j in np.arange( 0, L, L12 ):
                coords.append([ i,j ])

    elif ( dimensions == 3 ):
        # Get important properties
        L13 = L / NParticles ** (1/3)
        phi = NParticles / L ** 3 # Particle density

        for i in np.arange( 0, L, L13 ):
            for j in np.arange( 0, L, L13 ):
                for k in np.arange( 0, L, L13 ):
                    coords.append([ i,j,k ])

    coords = np.array(coords)

    print ("The particle density is %5.4f" % (phi) )
    np.savetxt("geometry0.xyz", coords)


    print ("The length of coordinate array is:", len(coords))
    return coords


if ( __name__ == "__main__" ):
    initParticles()
