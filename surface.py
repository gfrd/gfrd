import math
import numpy
import random

from shape import *
from single import *
from pair import *
from utils import *


class Surface( object ):
    """Surface should be added to the egfrd simulator by calling 
    sim.addSurface().

    """

    def __init__( self, name ):
        self.name = name


    def distanceTo( self, pos ):
        """Compute the absolute distance between pos and a closest point on 
        the surface.

        signedDistanceTo on the other hand returns a positive or negative 
        value depending on in which side the position is.  When this surface 
        defines a closed region in space, negative value means that the 
        position is inside.

        """

        return abs( self.signedDistanceTo( pos ) )
   

    def __str__( self ):
        return self.name


class PlanarSurface( Surface, Box ):
    """For example a membrane.

    Movement in 2D.

    """

    def __init__( self, origin, vectorX, vectorY, Lx, Ly, Lz=0, 
                  name="PlanarSurface" ):
        """Constructor.
        
        !origin! -- [ x0, y0, z0 ] is the center of the planar surface.
        vectorX -- [ x1, y1, z1 ] and
        vectorY -- [ x2, y2, z2 ] are 2 perpendicular vectors that don't have 
        to be normalized that span the plane. For example [1,0,0] and [0,1,0]
        for a plane at z=0.

        !Lx! -- lx and 
        !Ly! -- ly are the distances from the origin of the plane along vectorX 
            or vectorY to an edge of the plane. PlanarSurfaces are finite.
        Lz = dz, the thickness of the planar surface, can be omitted for Lz=0.

        """
        Surface.__init__( self, name )

        assert numpy.dot( vectorX, vectorY ) == 0.0
        # Orientation of surface is decided here.
        vectorZ = numpy.cross( vectorX, vectorY )
        Box.__init__( self, origin, vectorX, vectorY, vectorZ, Lx, Ly, Lz ) 
        self.defaultSingle = PlanarSurfaceSingle
        self.defaultPair = PlanarSurfacePair
        self.defaultInteractionSingle = PlanarSurfaceInteraction


    def drawBDdisplacement( self, dt, D ):
        r = math.sqrt( 2.0 * D * dt )
        # Draw 2 numbers from normal distribution.
        x, y = numpy.random.normal( 0.0, r, 2 )
        return x * self.unitX + y * self.unitY


    def randomVector( self, r ):
        x, y = randomVector2D( r )
        return x * self.unitX + y * self.unitY


    def randomPosition( self ):
        """Only uniform if vectorX and vectorY have same length.

        """
        return self.origin + random.uniform( -1, 1 ) * self.vectorX + \
                             random.uniform( -1, 1 ) * self.vectorY


    def minimalDistanceFromSurface( self, radius ):
        """A particle that is not on this surface has to be at least this far 
        away from the surface (measured from the origin of particle to the z = 
        0 plane of the surface).

        """
        return (self.Lz + radius) * MINIMAL_SEPERATION_FACTOR


    def randomUnbindingSite( self, pos, radius ):
        # Todo. SAFETY.
        return pos + random.choice( [ -1, 1 ] ) * \
                     self.minimalDistanceFromSurface( radius )  * self.unitZ


class CylindricalSurface( Surface, Cylinder ):
    """For example the DNA.

    Movement in 1D.

    """

    def __init__( self, origin, radius, orientation, size, 
                  name="CylindricalSurface" ):
        """Constructor.
        
        !origin! -- [ x0, y0, z0 ] is the center of the cylindrical surface.
        radius -- r is the radis of the cylinder.
        orientation -- [ x1, y1, z1 ] is a vector that doesn't have to
            normalized that defines the orienation of the cylinder. For 
            example [0,0,1] for a for a cylinder along the z-axis.
        !size! -- lz is the distances from the origin of the cylinder along 
            the oriention vector to the end of the cylinder. So effectively
            the half-length. CylindricalSurfaces are finite.

        """
        Surface.__init__( self, name )
        Cylinder.__init__( self, origin, radius, orientation, size )
        self.defaultSingle = CylindricalSurfaceSingle
        self.defaultPair = CylindricalSurfacePair
        self.defaultInteractionSingle = CylindricalSurfaceInteraction


    def drawBDdisplacement( self, dt, D ):
        r = math.sqrt( 2.0 * D * dt )
        # Draw 1 number from normal distribution.
        z = numpy.random.normal( 0.0, r, 1 )
        return z * self.unitZ


    def randomVector( self, r ):
        return random.choice( [-1, 1] ) * r * self.unitZ


    def randomPosition( self ):
        return self.origin + random.uniform( -1, 1 ) * self.vectorZ


    def minimalDistanceFromSurface( self, radius ):
        """A particle that is not on this surface has to be at least this far 
        away from the surface (measured from the origin of the particle to the 
        the central axis of the surface.

        """
        return ( self.radius + radius ) * MINIMAL_SEPERATION_FACTOR


    def randomUnbindingSite( self, pos, radius ):
        # Todo. SAFETY.
        x, y = randomVector2D( self.minimalDistanceFromSurface( radius ) )
        return pos + x * self.unitX + y * self.unitY


class CuboidalRegion( Surface, Box ):
    """
    A region that is (and can be) used for throwing in particles.

    Do not try to add this as a surface to your simulator, it won't work.

    It is also a Surface because particles for which no surface is 
    specified are tagged surface = defaultSurface, which is an instance of 
    this class. See gfrdbase.py.

    Movement in 3D.

    """

    def __init__( self, corner, size, name='world' ):
        """ Constructor.

        corner -- [ x0, y0, z0 ] is one corner of the cube.
        size -- [ sx, sy, sz ] is the vector from the origin to the diagonal
        point.

        """
        Surface.__init__( self, name )
        self.size = numpy.array( size )
        Lx = size[0]
        Ly = size[1]
        Lz = size[2]
        Box.__init__( self, corner + self.size / 2, [ Lx, 0, 0 ], [ 0, Ly, 0 ],
                      [ 0, 0, Lz ], Lx / 2, Ly / 2, Lz / 2 ) 
        self.defaultSingle = SphericalSingle
        self.defaultPair = SphericalPair


    def drawBDdisplacement( self, dt, D ):
        r = math.sqrt( 2.0 * D * dt )
        # Draw 3 numbers from normal distribution.
        return numpy.random.normal( 0.0, r, 3 )


    def randomVector( self, length ):
        return randomVector( length )


    def signedDistanceTo( self, pos ):
        """Overrule signedDistanceTo from Box. 
        Only for CuboidalRegions is cyclicTranspose 'pos' not needed.

        """
        raise RuntimeError( 'This method should not be used. Did you '
                            'accidently add this CuboidalRegion to the '
                            'surfacelist using s.addSurface()?' )
        corner = self.origin - size / 2
        dists = numpy.concatenate( ( corner - pos,
                                     corner + self.size - pos ) )
        absdists = numpy.abs( dists )
        i = numpy.argmin( absdists )
        return dists[i]


    def randomPosition( self ):
        """Returns a random position equidistributed within
        the region in space defined by negative signed distance.
        See also signedDistanceTo().

        """
        corner = self.origin - self.size / 2
        return numpy.array( [ random.uniform( corner[0], self.size[0] ),
                              random.uniform( corner[1], self.size[1] ),
                              random.uniform( corner[2], self.size[2] ) ] )


