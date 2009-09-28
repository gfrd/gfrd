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

    def __init__( self, origin, vectorX, vectorY, Lx, Ly, Lz=None, 
                  name="PlanarSurface" ):
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


    def minimalOffset( self, radius ):
        """A particle that is not on this surface has to be at least this far 
        away from the z = 0-plane of the surface.

        """
        return (self.Lz + radius) * UNBIND_SAFETY


    def randomUnbindingSite( self, pos, radius ):
        # Todo. SAFETY.
        return pos + random.choice( [ -1, 1 ] ) * \
                     self.minimalOffset( radius )  * self.unitZ


class CylindricalSurface( Surface, Cylinder ):
    """For example the DNA.

    Movement in 1D.

    """

    def __init__( self, origin, radius, orientation, size, 
                  name="CylindricalSurface" ):
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


    def minimalOffset( self, radius ):
        """A particle that is not on this surface has to be at least this far 
        away from the central axis of the surface.

        """
        return ( self.radius + radius ) * UNBIND_SAFETY


    def randomUnbindingSite( self, pos, radius ):
        # Todo. SAFETY.
        x, y = randomVector2D( self.minimalOffset( radius ) )
        return pos + x * self.unitX + y * self.unitY


class CuboidalSurface( Surface, Box ):
    """Surface that is only used for throwing in particles. Those particles 
    will than later be tagged with surface = defaultSurface, which is an 
    instance of the World class. See gfrdbase.py.

    If no surface is specified, particles are tagged with an instance of this 
    one.

    Movement in 3D.

    """

    def __init__( self, origin, size, name='world' ):
        """origin -- = [ x0, y0, z0 ] is one edge of the cube.
        size -- = [ sx, sy, sz ] is the vector from the origin to the diagonal
        point.

        """
        Surface.__init__( self, name )
        self.size = numpy.array( size )
        Lx = size[0]
        Ly = size[1]
        Lz = size[2]
        Box.__init__( self, origin + self.size / 2, [ Lx, 0, 0 ], [ 0, Ly, 0 ],
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
        Only for CuboidalSurfaces is cyclicTranspose 'pos' not needed.

        """
        raise RuntimeError( 'This method should not be used. Did you '
                            'accidently add this CuboidalSurface to the '
                            'surfacelist using s.addSurface()?' )
        edge = self.origin - size / 2
        dists = numpy.concatenate( ( edge - pos,
                                     edge + self.size - pos ) )
        absdists = numpy.abs( dists )
        i = numpy.argmin( absdists )
        return dists[i]


    def randomPosition( self ):
        """Returns a random position equidistributed within
        the region in space defined by negative signed distance.
        See also signedDistanceTo().

        """
        edge = self.origin - self.size / 2
        return numpy.array( [ random.uniform( edge[0], self.size[0] ),
                              random.uniform( edge[1], self.size[1] ),
                              random.uniform( edge[2], self.size[2] ) ] )


