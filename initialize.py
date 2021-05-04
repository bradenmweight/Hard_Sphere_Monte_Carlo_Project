import numpy as np

import Parameters

def getDists(coords):
    distx = np.subtract.outer(coords[:,1],coords[:,1])
    disty = np.subtract.outer(coords[:,2],coords[:,2])
    distz = np.subtract.outer(coords[:,3],coords[:,3])

    distr = np.sqrt( distx**2 + disty**2 + distz**2 )
    np.savetxt("Initial_Dists.dat", distr)

    print ("Min Distance:", distr[0,1] )
    print ("Box Size:", Parameters.Parameters.latticeLength)
    
def initParticles_old(  ):
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
    
    # Simple Cubic Lattice
    assert ( latticeType == "SC" and dimensions != 1 ), "Only 2D and 3D SC lattice is currently implemented."
    
    L = latticeLength
    radius = Parameters.Parameters.particleDiameter / 2
    V = L ** dimensions

    coords = []
    if ( dimensions == 2 ):
        L12 = L / NParticles ** (1/2)
        for i in np.arange( 0, L, L12 ):
            for j in np.arange( 0, L, L12 ):
                coords.append([ 1,i,j ])

    elif ( dimensions == 3 ):
        L13 = L / NParticles ** (1/3)
        sites = np.arange( 0, L, L13 )
        maxR = np.max( sites )
        for i in sites:
            for j in sites:
                for k in sites:
                    if ( maxR in [i,j,k] ):
                        continue
                    coords.append([ 1,i,j,k ])

    coords = np.array(coords)

    np.savetxt("geometry0.xyz", coords)

    
    if ( dimensions == 3 ):
        Vs = (4/3) * np.pi * radius**3
    elif ( dimensions == 2):
        Vs = np.pi * radius**2

    n =  len(coords) / V # Particle density
    phi = len(coords) * Vs / V # Volume Fraction

    print ("\n\tNumber Density:", n)
    print ("\n\tVolume Fraction:", phi)

    file01 = open("Density_data.dat","w")
    file01.write("Number_Density n (V^-1)\tVolume_Fraction \phi\n")
    file01.write("%5.5f\t%5.5f\n" % (n, phi))

    #print ("The length of coordinate array is:", len(coords))
    # Define boundary points
    #Lmin = np.array([ np.min(coords[:,d+1]) for d in range(dimensions) ]) #- Parameters.Parameters.particleDiameter
    #Lmax = np.array([ np.max(coords[:,d+1]) for d in range(dimensions) ]) #+ Parameters.Parameters.particleDiameter
    Lmin = np.array([0,0,0])
    Lmax = np.array([L,L,L])
    print ("Max and Min")
    print (Lmax)
    print (Lmin)

    #getDists(coords)




    return coords, Lmin, Lmax 

def initParticles_fromPHI0():
    """
    Initialize particles according to specified lattice arrangement and size.

    OUTPUT:
    coords:         nd.array() -- shape = NParticles x 4 = Type, x,y,z

    """
    
    #L = Parameters.Parameters.latticeLength
    #V = L ** dimensions

    latticeType = Parameters.Parameters.latticeType
    dimensions = Parameters.Parameters.dimensions
    radius = Parameters.Parameters.particleDiameter / 2
    phi0 = Parameters.Parameters.phi0
    NParticles = Parameters.Parameters.NParticles

    if ( dimensions == 3 ):
        Vs = (4/3) * np.pi * radius**3
        Ns = np.ceil(NParticles ** (1/3))
    elif ( dimensions == 2):
        Vs = np.pi * radius**2
        Ns = np.ceil(NParticles ** (1/2))

    # phi = len(coords) * Vs / V
    V = NParticles * Vs / phi0
    L = np.round( V ** (1/3) ,5)

    coords = []
    if ( dimensions == 2 ):
        L12 = L / (Ns**2+1) ** (1/2)
        sites = np.arange( 0, L, L12 )
        maxR = np.max(sites)
        for i in sites:
            for j in sites:
                if ( maxR in [i,j] ):
                        continue
                coords.append([ 1,i,j ])

    elif ( dimensions == 3 ):
        L13 = L / (Ns**3+1) ** (1/3)
        if ( latticeType == "SC" ):
            sites = np.arange( 0, L, L13 )
            if ( len(sites)**3 == NParticles ):
                print ("NEED MORE SITES!")
                exit()

            maxR = np.max(sites)
            print ("Total Sites:",len(sites)**3)
            print ("(N,Ns):",NParticles,Ns)
            for i in sites:
                for j in sites:
                    for k in sites:
                        if ( maxR in [i,j,k] or len(coords) >= NParticles ):
                            continue
                        coords.append([ 1,i,j,k ])

        elif ( latticeType == "FCC" ):
            sites = np.arange( 0, L, L13 )
            if ( len(sites)**3 == NParticles ):
                print ("NEED MORE SITES!")
                exit()
            maxR = np.max(sites)
            print ("Total Sites:",len(sites)**3)
            print ("(N,Ns):",NParticles,Ns)
            for i in sites:
                for j in sites:
                    for k in sites:
                        if ( maxR in [i,j,k] or len(coords) >= NParticles ):
                            continue
                        coords.append([ 1,i,j,k ])


 




    coords = np.array(coords)

    np.savetxt("geometry0.xyz", coords)

    
    n =  len(coords) / V # Particle density
    phi = len(coords) * Vs / V # Volume Fraction

    print ("\n\tNumber Density:", n)
    print ("\n\tVolume Fraction:", phi, phi0)

    file01 = open("Density_data.dat","w")
    file01.write("Number_Density n (V^-1)\tVolume_Fraction \phi\n")
    file01.write("%5.5f\t%5.5f\n" % (n, phi))

    #print ("The length of coordinate array is:", len(coords))
    # Define boundary points
    #Lmin = np.array([ np.min(coords[:,d+1]) for d in range(dimensions) ]) #- Parameters.Parameters.particleDiameter
    #Lmax = np.array([ np.max(coords[:,d+1]) for d in range(dimensions) ]) #+ Parameters.Parameters.particleDiameter
    Lmin = np.array([0,0,0])
    Lmax = np.array([L,L,L])
    print ("Max and Min")
    print (Lmax)
    print (Lmin)

    #getDists(coords)


    return coords, Lmin, Lmax 


if ( __name__ == "__main__" ):
    initParticles()

