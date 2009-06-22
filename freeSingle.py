from single import Single
from shape import Sphere
from domain import SimpleDomain
from utils import randomVector, length
from _gfrd import EventType, FirstPassageGreensFunction

# Spheres, Circles, Segments in free space.
# Maybe remove this layer of abstraction later.
class FreeSingle( Single ):
    def __init__( self, particle, reactionTypes, gf ):
        Single.__init__( self, particle, reactionTypes, \
                [ SimpleDomain( 0, 0, (None, 0), gf ) ] )

    def reset( self ):
        self.setSize( self.getMinSize() )
        self.dt = 0.0
        self.eventType = EventType.ESCAPE

    def isReset( self ):
        return self.size == self.getMinSize() and self.dt == 0.0\
               and self.eventType == EventType.ESCAPE

    def getPos( self ):
        return self.shellList[0].pos
    def setPos( self, pos ):
        self.shellList[0].pos = pos
        self.particle.pos = pos
    pos = property( getPos, setPos )

    def setSize( self, size ):
        assert size - self.getMinSize() >= 0.0
        self.shellList[0].size = size
        self.domains[0].a = self.getMobilitySize()
    def getSize( self ):
        return self.shellList[0].size
    size = property( getSize, setSize )

    def getMinSize( self ):
        return self.particle.species.radius

    def getMobilitySize( self ):
        return self.size - self.getMinSize()


class SphericalSingle( FreeSingle ):
    def __init__( self, particle, reactionTypes ):
        gf = FirstPassageGreensFunction( particle.species.D )
        FreeSingle.__init__( self, particle, reactionTypes, gf )
        self.shellList = [ Sphere( particle.pos, self.getMinSize() ), ]

    def toExternal( self, domains ):
        # Note: you need to be apply self.applyBoundary after calling this.
        r = domains[0]
        displacement = randomVector(r)
        assert abs( length( displacement ) - r ) <= 1e-15 * r
        newpos = self.pos + displacement
        return newpos


class CylindricalSingle1D( FreeSingle ):
    def __init__( self, particle, reactionTypes, cylindricalSurface ):
        radius = cylindricalSurface.outside.radius
        orientation = cylindricalSurface.outside.orientation
        gf = FirstPassageGreensFunction( particle.species.D )
        FreeSingle.__init__( self, particle, reactionTypes, gf )
        self.shellList = [ Cylinder( particle.pos, radius, orientation, self.getMinSize() ), ]

    def toExternal( self, domains ):
        # Note: you need to be apply self.applyBoundary after calling this.
        r = domains[0]
        displacement = r*orientation
        assert abs( length( displacement ) - r ) <= 1e-15 * r
        newpos = self.pos + displacement
        return newpos













################### Todo ##############################
class CircularSingle( FreeSingle ):
    def __init__( self ):
        pass


