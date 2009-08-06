import math
import numpy
from _gfrd import EventType, FirstPassageGreensFunction, FirstPassagePairGreensFunction
from utils import INF, NOWHERE, SAFETY, randomVector, randomVector2D, length, normalize, rotateVector
from shape import Sphere, Cylinder
from domain import RadialDomain1D, CartesianDomain1D, RadialDomain2D


'''
Not true anymore:
Note: a single is unaware of surface it is on.  Retrieve from particle if 
you need to know. (No duplication of information).
'''
class Single( object ):
    def __init__( self, particle, reactionTypes ):
        self.particle = particle
        self.surface = particle.surface
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


    def getPos( self ):
        return self.shellList[0].origin
    def setPos( self, pos ):
        self.shellList[0].origin = pos
        self.particle.pos = pos
    pos = property( getPos, setPos )


    def posString( self ):
        return '(%3.1f %3.1f %3.1f)' % ( self.pos[0], self.pos[1], self.pos[2] ) 


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
        return 'Single' + str( self.particle ) + '. pos=' + self.posString()


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


    def isReset( self ):
        return self.radius == self.getMinRadius() and self.dt == 0.0\
               and self.eventType == EventType.ESCAPE


    def __str__( self ):
        return 'Free' + Single.__str__( self )


class SphericalSingle3D( FreeSingle ):
    def __init__( self, particle, reactionTypes ):
        FreeSingle.__init__( self, particle, reactionTypes )

        self.shellList = [ Sphere( particle.pos, self.getMinRadius() ) ]

        # Create a radial domain of size mobilityRadius=0.
        gf = FirstPassageGreensFunction( particle.species.D )
        self.domains = [ RadialDomain1D( self.getMobilityRadius(), gf ) ]

    def calculateDisplacement( self, r ):
        return randomVector(r)


'''
Hockey pucks.
'''
class CylindricalSingle2D( FreeSingle ):
    def __init__( self, particle, reactionTypes ):
        FreeSingle.__init__( self, particle, reactionTypes )

        self.shellList = [ Cylinder( particle.pos, self.getMinRadius(), self.surface.unitZ, self.surface.Lz * SAFETY + 2 * particle.radius ) ]

        # Create a radial domain of size mobilityRadius=0.
        gf = FirstPassageGreensFunction( particle.species.D )
        self.domains = [ RadialDomain1D( self.getMobilityRadius(), gf ) ]


    def calculateDisplacement( self, r ):
        x, y = randomVector2D( r )
        return x * self.surface.unitX + y * self.surface.unitY


'''
Rods.
'''
class CylindricalSingle1D( FreeSingle ):
    def __init__( self, particle, reactionTypes ):
        FreeSingle.__init__( self, particle, reactionTypes )

        # Heads up. Cylinder's size is determined by getMinRadius().
        self.shellList = [ Cylinder( particle.pos, self.surface.radius * SAFETY + 2 * particle.radius, self.surface.unitZ, self.getMinRadius() ) ]

        # Create a cartesian domain of size mobilityRadius=0.
        gf = FirstPassageGreensFunction( particle.species.D )
        self.domains = [ CartesianDomain1D( 0, (0, 0), self.getMobilityRadius(), gf ) ]


    def calculateDisplacement( self, z ):
        return z * self.shellList[0].unitZ


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


'''
Singles close to a surface.
'''
class InteractionSingle( Single ):
    def __init__( self, particle, surface, reactionTypes, interactionType ):
        self.interactionSurface = surface
        self.interactionType = interactionType
        Single.__init__( self, particle, reactionTypes )


    def reset( self ):
        pass


    def __str__( self ):
        return 'Interaction' + Single.__str__( self )



'''
Interaction with plane.
'''
class InteractionSingle2D( InteractionSingle ):
    def __init__( self, particle, surface, reactionTypes, interactionType, origin, radius, orientationVector, size, particleOffset, projectedPoint = None ):
        InteractionSingle.__init__( self, particle, surface, reactionTypes, interactionType )

        self.shellList = [ Cylinder( origin, radius, orientationVector, size ) ]

        # Free diffusion in r direction.
        gfr = FirstPassageGreensFunction( particle.species.D )
        rDomain = RadialDomain1D( radius - particle.species.radius, gfr )

        # Interaction possible in z direction.
        # Todo. Correct gf.
        gfz = FirstPassageGreensFunction( particle.species.D )
        zDomain = CartesianDomain1D( particleOffset[1], (interactionType.k, 0), size - particle.species.radius, gfz )

        self.domains = [ rDomain, zDomain ]


    def drawNewPosition( self, dt ):
        r = self.domains[0].drawPosition( dt )
        z = self.domains[1].drawPosition( dt )
        x, y = randomVector2D( r )
        return self.pos + x * self.interactionSurface.unitX + y * self.interactionSurface.unitY + z * self.shellList[0].unitZ


'''
Interaction with cylinder.
'''
class InteractionSingle1D( InteractionSingle ):
    def __init__( self, particle, surface, reactionTypes, interactionType, origin, radius, orientationVector, size, particleOffset, projectedPoint ):
        self.unitR = normalize( particle.pos - projectedPoint ) # Not needed for InteractionSingle2D
        InteractionSingle.__init__( self, particle, surface, reactionTypes, interactionType )

        self.shellList = [ Cylinder( origin, radius, orientationVector, size ) ]

        # Interaction possible in r direction.
        # Todo. Correct gf.
        gfr = FirstPassagePairGreensFunction( particle.species.D, interactionType.k, surface.radius )
        self.pgf = gfr
        gfr.seta( radius )
        rDomain = RadialDomain2D( surface.radius + particle.species.radius, particleOffset[0], radius - particle.species.radius, gfr )

        # Free diffusion in z direction.
        gfz = FirstPassageGreensFunction( particle.species.D )
        zDomain = CartesianDomain1D( particleOffset[1], (0, 0), size - particle.species.radius, gfz )

        self.domains = [ rDomain, zDomain ]


    def drawNewPosition( self, dt ):
        # Todo, gf.
        gf = self.choosePairGreensFunction( dt )
        r, theta = self.domains[0].drawPosition( gf, dt )
        z = self.domains[1].drawPosition( dt )
        # Calculate new position starting from origin.
        newVectorR = r * rotateVector( self.unitR, self.shellList[0].unitZ, theta )
        return self.pos + newVectorR + z * self.shellList[0].unitZ


    def choosePairGreensFunction( self, dt ):
        # Todo.
        return self.pgf

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

