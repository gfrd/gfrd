import math
import numpy
from _gfrd import EventType, FirstPassageGreensFunction
from utils import INF, NOWHERE, randomVector, randomVector2D, length
from shape import Sphere, Cylinder
from domain import RadialDomain1D, CartesianDomain1D


'''
Note: a single is unaware of surface it is on.  Retrieve from particle if 
you need to know. (No duplication of information).
'''
class Single( object ):
    def __init__( self, particle, reactionTypes ):
        self.particle = particle
        self.reactionTypes = reactionTypes
        if reactionTypes:
            self.k_tot = sum( ( rt.k for rt in reactionTypes ) )
        else: 
            self.k_tot = 0
        self.eventType = None
        self.eventID = None
        self.multiplicity = 1


    def getD( self ):
        return self.particle.species.D


    def initialize( self, t ):
        self.reset()
        self.lastTime = t
        self.domains


    '''
    Returns an (escapeTime, eventType, activeDomain)-tuple.
    By returning the arguments it is a pure function. 
    '''
    def determineNextEvent( self ):
        return min(self.drawEscapeOrInteractionTime(), 
                self.drawReactionTime()) 


    '''
    Returns an (escapeTime, eventType, activeDomain)-tuple.
    Handles also all interaction events.
    '''
    def drawEscapeOrInteractionTime( self ):
        if self.getD() == 0:
            return INF, EventType.ESCAPE, None
        else:
            '''
            Note: we are not calling domain.drawEventType() just yet, but 
            postpone it to the very last minute (when this event is executed 
            in fireSingle), and memorize the activeDomain like this.

            So this can still be an interaction or an escape.

            Also note that in case this single will get a reaction event 
            instead of this escape event (its dt is smaller in 
            determineNextEvent), self.activeDomain won't be used at all, as 
            long as you make sure reaction events are taken care of before 
            escape events in fireSingle.

            By not letting the domains not notice that one of them has been 
            made active, we also don't have to reset a flag or something when
            there is a burst (again, just ignore the activeDomain flag).
            '''
            return min( (d.drawTime(), EventType.ESCAPE, d) for d in self.domains )


    '''
    Returns a (reactionTime, eventType, activeDomain=None)-tuple.
    '''
    def drawReactionTime( self ):
        if self.k_tot == 0:
            dt = INF
        elif self.k_tot == INF:
            dt = 0.0
        else:
            dt = ( 1.0 / self.k_tot ) * math.log( 1.0 / numpy.random.uniform() )
        return dt, EventType.REACTION, None


    '''
    Copy pasted.
    '''
    def drawReactionType( self ):
        k_array = [ rt.k for rt in self.reactionTypes ]
        k_array = numpy.add.accumulate( k_array )
        k_max = k_array[-1]

        rnd = numpy.random.uniform()
        i = numpy.searchsorted( k_array, rnd * k_max )
        return self.reactionTypes[i]


    def check( self ):
        pass


    def __str__( self ):
        return 'Shell' + str( self.particle ) + (" (%3.1f %3.1f %3.1f)" % (self.pos[0], self.pos[1], self.pos[2]) )


'''
Singles in free space.
Maybe remove this layer of abstraction later.
'''
class FreeSingle( Single ):
    def __init__( self, particle, reactionTypes ):
        Single.__init__( self, particle, reactionTypes )


    def drawNewPosition( self, dt ):
        r = self.domains[0].drawPosition( dt )
        displacement = self.calculateDisplacement( r )
        assert abs( length( displacement ) - abs(r) ) <= 1e-15 * abs(r)
        return self.pos + displacement


    def reset( self ):
        self.radius = self.getMinRadius()
        self.dt = 0.0
        self.eventType = EventType.ESCAPE
      
        # Todo. Cant't we set proper size and event immediately? (like with 
        # pairs).
        # Just set some activeDomain, doesn't matter which.
        self.activeDomain = self.domains[0] # FIXME


    def isReset( self ):
        return self.radius == self.getMinRadius() and self.dt == 0.0\
               and self.eventType == EventType.ESCAPE


    def getPos( self ):
        return self.shellList[0].origin
    def setPos( self, pos ):
        self.shellList[0].origin = pos
        self.particle.pos = pos
    pos = property( getPos, setPos )


    def getRadius( self ):
        return self.shellList[0].radius
    def setRadius( self, radius ):
        assert radius - self.getMinRadius() >= 0.0
        self.shellList[0].radius = radius
        '''
        A bit tricky: getMobilityRadius() uses self.radius, which is 
        getRadius(), which is self.shellList[0], which is already updated.
        Works for cylindricalSingle1D as well!
        '''
        a = self.getMobilityRadius()
        self.domains[0].a = a
    radius = property( getRadius, setRadius )


    def getMinRadius( self ):
        # Also works for cylinders (2D/3D).
        return self.particle.species.radius


    def getMobilityRadius( self ):
        return self.radius - self.getMinRadius()


class SphericalSingle3D( FreeSingle ):
    def __init__( self, particle, reactionTypes, distFunc ):
        FreeSingle.__init__( self, particle, reactionTypes )

        self.shellList = [ Sphere( particle.pos, self.getMinRadius(), distFunc ) ]

        # Create a radial domain of size mobilityRadius=0.
        gf = FirstPassageGreensFunction( particle.species.D )
        self.domains = [ RadialDomain1D( self.getMobilityRadius(), gf ) ]

    def calculateDisplacement( self, r ):
        return randomVector(r)


class CylindricalSingle2D( FreeSingle ):
    def __init__( self, particle, reactionTypes, distFunc ):
        FreeSingle.__init__( self, particle, reactionTypes )

        self.shellList = [ Cylinder( particle.pos, self.getMinRadius(), particle.surface.normal, particle.radius, distFunc ) ]

        # Create a radial domain of size mobilityRadius=0.
        gf = FirstPassageGreensFunction( particle.species.D )
        self.domains = [ RadialDomain1D( self.getMobilityRadius(), gf ) ]


    def calculateDisplacement( self, r ):
        x, y = randomVector2D( r )
        return x * self.particle.surface.outside.xUnitVector + y * self.particle.surface.outside.yUnitVector


class CylindricalSingle1D( FreeSingle ):
    def __init__( self, particle, reactionTypes, distFunc ):
        FreeSingle.__init__( self, particle, reactionTypes )

        # Heads up. Cylinder's size is determined by getMinRadius().
        self.shellList = [ Cylinder( particle.pos, particle.radius, particle.surface.outside.orientation, self.getMinRadius(), distFunc ) ]

        # Create a cartesian domain of size mobilityRadius=0.
        gf = FirstPassageGreensFunction( particle.species.D )
        self.domains = [ CartesianDomain1D( 0, (0, 0), self.getMobilityRadius(), gf ) ]


    def calculateDisplacement( self, r ):
        return r * self.shellList[0].orientation


    '''
    Overloaded methods getRadius and setRadius. Nice trick. property() needs 
    to be redefined as well.
    '''
    def getRadius( self ):
        # Heads up. Return cylinder's size.
        return self.shellList[0].size
    def setRadius( self, size ):
        assert size - self.getMinRadius() >= 0.0 # Still fine.
        # Heads up. A larger shell means a larger CylindricalSingle2D's size.
        self.shellList[0].size = size
        self.domains[0].a = self.getMobilityRadius()
    radius = property( getRadius, setRadius )


class DummySingle( object ):
    def __init__( self ):
        self.multiplicity = 1
        self.size = 0.0
        self.shellList = [ Sphere( NOWHERE, 0.0 ), ]


    def getMinSize( self ):
        return 0.0


    def getD( self ):
        return 0.0


    def getPos( self ):
        return NOWHERE
    pos = property( getPos )


    def __str__( self ):
        return 'DummySingle()'

