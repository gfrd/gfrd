import scipy
import math
import numpy

Pi = scipy.pi
Pi2 = scipy.pi * 2.0
PiSqrt = math.sqrt( scipy.pi )

N_A = 6.0221367e23
INF = numpy.inf

ZEROPOS = numpy.array( [ 0., 0., 0. ] )
NOWHERE = numpy.array( ( INF, INF, INF ) )


### Floating points
# Safetey used for not making sure shells do not touch each other.
# SAFETY should be smaller than MINIMAL_SINGLE_RADIUS_FACTOR.
SAFETY = 1.0 + 1e-5

# Tolerance used for float comparison functions. Oversimplifying: two floats a 
# and b are considered to be equal if abs( a - b ) < TOLERANCE * abs( a ).
TOLERANCE = 1e-8

# Multiplication factor used for seperating 2 particles or a particle and a 
# surface after unbinding.
MINIMAL_SEPERATION_FACTOR = 1.0 + TOLERANCE


### Surfaces
# Scan for surfaces in a spherical region with a radius of the particle radius 
# multiplied by INTERACTION_HORIZON_FACTOR.
# Should be bigger than the burst factors below, we check for interactions 
# with surfaces first.
# Should not be too big, because we can not handle 2 surfaces at once.
INTERACTION_HORIZON_FACTOR = 2

# Instead of letting spheres in the cytoplasm almost touch the surfaces, 
# mulitply the radii of the surfaces by SURFACE_SAFETY_ZONE to create a safety 
# zone for the particles on the surface.
# Todo. Not yet implemented.
SURFACE_SAFETY_ZONE = 1.1
''' Pseudocode from board. Todo.

shellsize = calculateSingleShellSize()
shellsize = min( shellsize, distanceToSurface )

if shellsize < MINIMAL_SINGLE_RADIUS_FACTOR * particleRadius
    # multi
    pass
else:
    # Reduce shellsize to allow for safety zone.
    shellsize -= SURFACE_SAFETY_ZONE * particleRadius
    shellsize = max( shellsize, minshellsize )
    # Make Single with radius 'shellsize'.
    pass
'''


### Bursting
# When forming a cylinder, clear a cylindrical region with a radius of the 
# particle radius multiplied by
# Todo. Not yet implemented.
CYLINDER_BURST_RADIUS_FACTOR = 1.1

# When forming a cylinder, clear a cylindrical region with a size of the 
# particle radius multiplied by
# Todo. Not yet implemented.
CYLINDER_BURST_SIZE_FACTOR = 1.1

# When forming a sphere, clear a spherical region with a radius of the 
# particle radius multiplied by
SPHERE_BURST_RADIUS_FACTOR = 1.1


### Minimal radii
# Multiplication factor that determines the minimal shell around a particle, 
# before falling back to brute force brownian dynamics. 
MINIMAL_SINGLE_RADIUS_FACTOR = 1.1

# Multiplication factor that determine the radius of the sphere around a 
# particle that is updated using brownian dynamics.
# Should be smaller than MINIMAL_SINGLE_RADIUS_FACTOR.
MULTI_SHELL_FACTOR = 1.05


### Temporary.
# Temporary constant to make sure initial position within a domain is not too 
# close to one of the boundaries.
MAX_DOMAIN_SIZE_FACTOR = 20



