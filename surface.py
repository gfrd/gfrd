import math
import numpy
import random

from shape import *
from single import *
from pair import *
from utils import *


'''
Surface should be added to the egfrd simulator by calling sim.addSurface().
'''
class Surface( object ):
    def __init__( self, name ):
        self.name = name


    '''
    This method computes the absolute distance between pos and
    a closest point on the surface.

    signedDistanceTo on the other hand returns a positive
    or negative value depending on in which side the position is.
    When this surface defines a closed region in space,
    negative value means that the position is inside.
    '''
    def distanceTo( self, pos ):
        return abs( self.signedDistanceTo( pos ) )
   

    def __str__( self ):
        return self.name


'''
For example a membrane.
'''
class PlanarSurface( Surface, Box ):
    def __init__( self, origin, vectorX, vectorY, Lx, Ly, Lz=None, name="PlanarSurface" ):
        Surface.__init__( self, name )

        assert numpy.dot( vectorX, vectorY ) == 0.0
        # Orientation of surface is decided here.
        vectorZ = numpy.cross( vectorX, vectorY )
        Box.__init__( self, origin, vectorX, vectorY, vectorZ, Lx, Ly, Lz ) 
        self.defaultSingle = PlanarSurfaceSingle
        self.defaultPair = PlanarSurfacePair
        self.defaultInteractionSingle = PlanarSurfaceInteraction


    def drawBDdisplacement( self, t, D ):
        ro = math.sqrt( 2.0 * D * t )
        x, y = numpy.random.normal( 0.0, ro, 2 )
        return x*self.unitX + y*self.unitY


    def randomVector( self, r ):
        x, y = randomVector2D( r )
        return x * self.unitX + y * self.unitY


    '''
    Only uniform if vectorX and vectorY have same length.
    '''
    def randomPosition( self ):
        return self.origin + random.uniform(-1,1)*self.vectorX + random.uniform(-1,1)*self.vectorY


    def randomUnbindingSite( self, pos, radius ):
        # Todo. SAFETY.
        return pos + random.choice( [-1,1] ) * (self.Lz + radius ) * SAFETY * self.unitZ


'''
For example the DNA.
'''
class CylindricalSurface( Surface, Cylinder ):
    def __init__( self, origin, radius, orientation, size, name="CylindricalSurface" ):
        Surface.__init__( self, name )
        Cylinder.__init__( self, origin, radius, orientation, size )
        self.defaultSingle = CylindricalSurfaceSingle
        self.defaultPair = CylindricalSurfacePair
        self.defaultInteractionSingle = CylindricalSurfaceInteraction


    def drawBDdisplacement( self, t, D ):
        ro = math.sqrt( 2.0 * D * t )
        z = numpy.random.normal( 0.0, ro, 1 )
        return z*self.unitZ


    def randomVector( self, r ):
        return random.choice( [-1, 1] ) * r * self.unitZ


    def randomPosition( self ):
        return self.origin + random.uniform(-1, 1) * self.vectorZ


    def randomUnbindingSite( self, pos, radius ):
        # Todo. SAFETY.
        x, y = randomVector2D( ( self.radius + radius ) * SAFETY )
        return pos + x * self.unitX + y * self.unitY


'''
Surface that is only used for throwing in particles. Those particles will than 
later be tagged with surface=defaultSurface, which is an instance of the World 
class. See gfrdbase.py.
'''
class CuboidalSurface( Surface, Box ):
    '''
    origin -- = [ x0, y0, z0 ] is one edge of the cube.
    size -- = [ sx, sy, sz ] is the vector from the origin to the diagonal
    point.
    '''
    def __init__( self, origin, size, name='CuboidalSurface' ):
        Surface.__init__( self, name )
        self.size = numpy.array(size)
        Lx = size[0]
        Ly = size[1]
        Lz = size[2]
        Box.__init__( self, origin + self.size/2, [Lx, 0, 0], [0, Ly, 0], [0, 0, Lz], Lx/2, Ly/2, Lz/2 ) 
 

    '''
    Overrule signedDistanceTo from Box. 
    Only for CuboidalSurfaces is cyclicTranspose 'pos' not needed.
    '''
    def signedDistanceTo( self, pos ):
        raise 'You should not add CuboidalSurfaces to the surfacelist using s.addSurface()!'
        edge = self.origin - size/2
        dists = numpy.concatenate( ( edge - pos,
                                     edge+self.size-pos ) )
        absdists = numpy.abs( dists )
        i = numpy.argmin( absdists )
        return dists[i]


    '''
    Returns a random position equidistributed within
    the region in space defined by negative signed distance.
    See also signedDistanceTo().
    '''
    def randomPosition( self ):
        edge = self.origin - self.size/2
        return numpy.array( [ random.uniform( edge[0], self.size[0] ),
                              random.uniform( edge[1], self.size[1] ),
                              random.uniform( edge[2], self.size[2] )])


'''
If no surface is specified, particles are tagged with an instance of this one.
'''
class World( Surface ):
    def __init__( self ):
        Surface.__init__( self, 'world' )
        self.defaultSingle = SphericalSingle
        self.defaultPair = SphericalPair


    def randomVector( self, length ):
        return randomVector( length )


    def drawBDdisplacement( self, t, D ):
        ro = math.sqrt( 2.0 * D * t )
        return numpy.random.normal( 0.0, ro, 3 )


