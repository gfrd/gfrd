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
# SAFETY should be smaller than MINIMAL_SINGLE_RADIUS_FACTOR. Todo.
SAFETY = 1.0 + 1e-5 #5e-2 is too big! Todo.

# Tolerance used for float comparison functions. Oversimplifying: two floats a 
# and b are considered to be equal if abs( a - b ) < TOLERANCE * abs( a ).
TOLERANCE = 1e-8

# Multiplication factor used for seperating 2 particles or a particle and a 
# surface after unbinding.
MINIMAL_SEPERATION_FACTOR = 1.0 + TOLERANCE


### Surfaces
# Multiplication factor for determining how close a particle should be to a 
# surface for an interaction to become feasible. Should not be too big, 
# because we can not handle 2 surfaces at once.
INTERACTION_HORIZON_FACTOR = 2


### Bursting
# Todo. When forming a cylinder, clear volume.
CYLINDER_BURST_RADIUS_FACTOR = 1.1

# Todo. When forming a cylinder, clear volume.
CYLINDER_BURST_SIZE_FACTOR = 1.1

# Multiplication factor for determining which shells to burst when trying to 
# form a pair.  
SPHERE_BURST_RADIUS_FACTOR = 1.1


### Minimal radii
# Multiplication factor that determines the minimal shell around a particle, 
# before falling back to brute force brownian dynamics. 
MINIMAL_SINGLE_RADIUS_FACTOR = 1.1

# Multiplication factor that determine the radius of the sphere around a 
# particle that is updated using brownian dynamics.
# Should be smaller than MINIMAL_SINGLE_RADIUS_FACTOR? Todo.
MULTI_SHELL_FACTOR = 1.05


### Todo
# Temporary constant to make sure initial position within a domain is not too 
# close to one of the boundaries.
MAX_DOMAIN_SIZE_FACTOR = 20

# Todo. Instead of letting spheres in the cytoplasm almost touch the surfaces, 
# mulitply by the radii of the surfaces by SURFACE_SAFETY_ZONE to create a 
# safety zone for the particles on the surface.
SURFACE_SAFETY_ZONE = 1.1


