from single import Single
from shape import Sphere, Cylinder
from domain import SimpleDomain
from utils import randomVector, length
from _gfrd import EventType, FirstPassageGreensFunction

# Spheres, Circles, Segments in free space.
# Maybe remove this layer of abstraction later.
class FreeSingle( Single ):
    def __init__( self, particle, shell, reactionTypes, gf ):
        # Create a domain with mobilitySize 0.
        domain = [ SimpleDomain( 0, 0, (None, 0), gf ) ]
        Single.__init__( self, particle, shell, reactionTypes, domain )


    def reset( self ):
        self.setSize( self.getMinSize() )
        self.dt = 0.0
        self.eventType = EventType.ESCAPE


    def isReset( self ):
        return self.size == self.getMinSize() and self.dt == 0.0\
               and self.eventType == EventType.ESCAPE


    def getPos( self ):
        return self.shellList[0].origin
    def setPos( self, pos ):
        self.shellList[0].origin = pos
        self.particle.pos = pos
    pos = property( getPos, setPos )


    def setSize( self, size ):
        assert size - self.getMinSize() >= 0.0
        self.shellList[0].size = size
        # Works for cylindricalSingle1D as well!
        self.domains[0].size = self.getMobilitySize()
    def getSize( self ):
        return self.shellList[0].size
    size = property( getSize, setSize )


    def getMinSize( self ):
        return self.particle.species.radius


    def getMobilitySize( self ):
        return self.size - self.getMinSize()



class SphericalSingle( FreeSingle ):
    def __init__( self, particle, reactionTypes, distFunc ):
        minSize = particle.radius # self.getMinSize()
        shell = Sphere( particle.pos, minSize, distFunc )
        gf = FirstPassageGreensFunction( particle.species.D )
        FreeSingle.__init__( self, particle, shell, reactionTypes, gf )


    def toExternal( self, domains ):
        # Note: you need to apply self.applyBoundary after calling this.
        r = domains[0]
        displacement = randomVector(r)
        assert abs( length( displacement ) - r ) <= 1e-15 * r
        newpos = self.pos + displacement
        return newpos



class CylindricalSingle1D( FreeSingle ):
    def __init__( self, particle, reactionTypes, distFunc ):
        assert particle.radius <= particle.surface.outside.radius # Particle is absorbed in the DNA for now.
        minSize = particle.radius # self.getMinSize()
        shell = Cylinder( particle.pos, particle.radius, particle.surface.outside.orientation, minSize, distFunc )
        gf = FirstPassageGreensFunction( particle.species.D )

        # Todo. Split up radial and cartesian Singles after all?
        domain = [ SimpleDomain( 0, 0, (0, 0), gf ) ]
        Single.__init__( self, particle, shell, reactionTypes, domain  )


    def toExternal( self, domains ):
        # Note: you need to apply self.applyBoundary after calling this.
        r = domains[0]
        displacement = r*self.shellList[0].orientation
        assert abs( length( displacement ) - abs(r) ) <= 1e-15 * abs(r)
        newpos = self.pos + displacement
        return newpos













################### Todo ##############################
class CircularSingle( FreeSingle ):
    def __init__( self ):
        pass


