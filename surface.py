#from utils import *
import math
import numpy
import random

from shape import Cylinder, Box
from freeSingle import SphericalSingle, CylindricalSingle1D, CylindricalSingle2D
from utils import length, normalize

class Surface( object ):
    def __init__( self, name ):
        self.name = name


    '''
    This method computes the absolute distance between pos and
    a closest point on the surface.
    '''
    def distanceTo( self, pos ):
        return abs( self.signedDistanceTo( pos ) )
   

    '''
    Same as distanceTo(), but the return value is positive
    or negative depending on in which side the position is.
    
    When this surface defines a closed region in space,
    negative value means that the position is inside.
    '''
    def signedDistanceTo( self, pos ):
        return self.outside.signedDistanceTo( pos )


    def __str__( self ):
        return self.name



class DefaultSurface( Surface ):
    # If no surface is specified, particles are tagged with this one.
    def __init__( self ):
        Surface.__init__( self, 'defaultSurface' )
        self.defaultSingle = SphericalSingle



class CylindricalSurface( Surface ):
    def __init__( self, origin, radius, orientation, size, name="CylindricalSurface" ):
        Surface.__init__( self, name )
        self.outside = Cylinder( origin, radius, orientation, size )
        self.defaultSingle = CylindricalSingle1D


    def randomPosition( self ):
        return self.outside.origin + random.uniform(-1, 1) * self.outside.orientation * self.outside.size



class PlanarSurface( Surface ):
    def __init__( self, origin, xVector, yVector, Lx, Ly, Lz=None, name="PlanarSurface" ):
        Surface.__init__( self, name )

        assert numpy.dot( xVector, yVector ) == 0 # Todo. To doubles.
        zVector = numpy.cross( xVector, yVector )
        self.outside = Box( origin, xVector, yVector, zVector, Lx, Ly, Lz ) 
        self.defaultSingle = CylindricalSingle2D

        # Hessian normal form [ nx, ny, nz, p ].
        self.normal = self.outside.zUnitVector # Extra.
        self.hessianNormal = [ self.normal[0], self.normal[1], self.normal[2], length(self.outside.origin) ]



    def signedDistance( self, pos ):
        return numpy.dot( self.hessianNormal,\
                  numpy.array( [ pos[0], pos[1], pos[2], 1.0 ] ) ) - self.outside.Lz


    def randomPosition( self ):
        return self.outside.origin + random.uniform(-1,1)*self.outside.xVector + random.uniform(-1,1)*self.outside.yVector



class CuboidalSurface( Surface ):
    '''
    origin -- = [ x0, y0, z0 ] is the origin.
    size -- = [ sx, sy, sz ] is the vector from the origin to the diagonal
              point.
    '''
    def __init__( self, origin, size, name='someCuboidalSurface' ):
        Surface.__init__( self, name )
        self.setParams( origin, size )
        self.defaultSingle = SphericalSingle
 

    def signedDistanceTo( self, pos ):
        dists = numpy.concatenate( ( self.origin - pos,
                                     self.origin+self.size-pos ) )
        absdists = numpy.abs( dists )
        i = numpy.argmin( absdists )

        return dists[i]


    '''
    Returns a random position equidistributed within
    the region in space defined by negative signed distance.
    See also signedDistance().
    '''
    def randomPosition( self ):
        return numpy.array( [ random.uniform( self.origin[0], self.size[0] ),
                              random.uniform( self.origin[1], self.size[1] ),
                              random.uniform( self.origin[2], self.size[2] ) ]
                            )


    def setParams( self, origin, size ):

        self.origin = numpy.array( origin )
        self.size = size
        self.midpoint = ( self.size - self.origin ) * 0.5

    def getParams( self ):
        return self.params

    

"""
class SphericalSurface( Surface ):

    '''
    origin -- = [ x0, y0, z0 ] is the origin and
    radius -- the radius of the sphere.

    ( x - x0 )^2 + ( y - y0 )^2 + ( z - z0 )^2 = r^2
    '''
    def __init__( self, origin, radius ):
        self.setParams( origin, radius )

    def signedDistanceTo( self, pos ):
        return math.sqrt( ( ( pos - self.origin ) ** 2 ).sum() ) - self.radius

    def randomPosition( self ):
        pos = randomUnitVectorS()
        pos[0] *= self.radius
        return sphericalToCartesian( pos )

    def setParams( self, origin, radius ):
        self.origin = numpy.array( origin )
        self.radius = radius

    def getParams( self ):
        return self.params





######## GRAVEYARD ##############

class OneDimensionalSurface( Surface ):

    '''
    Infinetely long line along z-axis.
    '''

    def __init__( self, origin, orientation = [ 0, 0, 1 ] ):
        self.setParams( origin, orientation )

    #def signedDistanceTo( self, pos ):
    #    return math.sqrt( (pos[0] - self.origin[0])** 2 + (pos[1] - self.origin[1])**2 )

    def randomPosition( self ):
        raise 'randomPosition not supported. specified space has infinite volume.'

    def setParams( self, origin, orientation ):
        self.origin = numpy.array( origin )
        self.orientation = orientation

    def getParams( self ):
        return self.params

    # Given a position in 3D coordinates, how far are we to the left or right
    # from the origin of the line.
    def position2location( self, pos ):

        posFromOrigin = numpy.subtract( pos, self.origin )
        return numpy.dot( self.orientation, posFromOrigin )
    
    # Given a location on the line, what is the position in 3D coordinates.
    def location2Position( self, location ):
        
        return numpy.multiply( location, self.orientation )
"""
