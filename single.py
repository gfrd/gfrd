import math
import numpy
from _gfrd import *
from utils import *
from shape import *
from domain import *


'''
There are 2 main types of Singles:
    * NonInteractionSingle
    * InteractionSingle (when the particle is nearby a surface)

Each type of Single defines a list of domains, see domain.py. For each domain the Green's function is specified.
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


    def posString( self ):
        factor = 1
        return '(%g %g %g)' % ( self.pos[0]*factor, self.pos[1]*factor, self.pos[2]*factor ) 


    def getMinRadius( self ):
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
            determineNextEvent), and even though activeDomain is set, it won't 
            be used at all, since reaction events are taken care of before 
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
        return str( self.particle ) + '. pos=' + self.posString()


##############################################################################
'''
1 Particle inside a shell, no other particles around. 

There are 3 types of NonInteractionSingles:
    * SphericalSingle: spherical shell, 3D movement.
    * PlanarSurfaceSingle: cylindrical shell, 2D movement.
    * CylindricalSurfaceSingle: cylindrical shell, 1D movement.
'''
class NonInteractionSingle( Single ):
    def __init__( self, particle, reactionTypes ):
        Single.__init__( self, particle, reactionTypes )


    def getPos( self ):
        return self.shellList[0].origin
    def setPos( self, pos ):
        # setPos is called when NonInteractionSingles are updated.
        self.shellList[0].origin = pos
        self.particle.pos = pos
    pos = property( getPos, setPos )


    def getRadius( self ):
        #raise RuntimeError, 'Todo'
        return self.shellList[0].radius
    def setRadius( self, radius ):
        #raise RuntimeError, 'Todo'
        # setRadius is called when NonInteractionSingles are updated.
        assert radius - self.getMinRadius() >= 0.0
        self.shellList[0].radius = radius
        # A bit tricky: getMobilityRadius() uses self.radius, which uses 
        # getRadius(), which uses self.shellList[0], which is already updated.
        self.domains[0].a = self.getMobilityRadius()
    radius = property( getRadius, setRadius )


    def drawNewPosition( self, dt ):
        r = self.domains[0].drawPosition( dt )
        displacement = self.displacement( r )
        assert abs( length( displacement ) - abs(r) ) <= 1e-15 * abs(r)
        return self.pos + displacement


    def initialize( self, t ):
        self.reset()
        self.lastTime = t


    def reset( self ):
        self.radius = self.getMinRadius()
        self.dt = 0.0
        self.eventType = EventType.ESCAPE


    def isReset( self ):
        return self.radius == self.getMinRadius() and self.dt == 0.0\
               and self.eventType == EventType.ESCAPE


'''
1 Particle inside a (spherical) shell not on any surface.

    * Particle coordinate inside shell: r,theta,phi.
    * Domain: radial r.
    * Initial position: r=0.
    * Selected randomly when drawing displacement vector: theta, phi.
'''
class SphericalSingle( NonInteractionSingle ):
    def __init__( self, particle, reactionTypes ):
        NonInteractionSingle.__init__( self, particle, reactionTypes )

        self.shellList = [ Sphere( particle.pos, self.getMinRadius() ) ]

        # Create a radial domain of size mobilityRadius=0.
        gf = FirstPassageGreensFunction( self.getD() )
        self.domains = [ RadialDomain( self.getMobilityRadius(), gf ) ]

    def displacement( self, r ):
        return randomVector(r)


    def __str__( self ):
        return 'SphericalSingle' + Single.__str__( self )



'''
1 Particle inside a (cylindrical) shell on a PlanarSurface. (Hockey pucks).

    * Particle coordinates on surface: x,y.
    * Domain: radial r. (determines x and y together with theta).
    * Initial position: r=0.
    * Selected randomly when drawing displacement vector: theta.
'''
class PlanarSurfaceSingle( NonInteractionSingle ):
    def __init__( self, particle, reactionTypes ):
        NonInteractionSingle.__init__( self, particle, reactionTypes )

        '''
        Hockey pucks stick a bit out of the surface, so that if the 
        biggest particle would undergo an unbinding reaction, it can still be 
        placed within the hockey puck (and thus it won't interfere with other 
        shells.
        '''
        self.shellList = [ Cylinder( particle.pos, self.getMinRadius(), self.surface.unitZ, self.surface.Lz * SAFETY + 2 * particle.radius ) ]

        # Create a radial domain of size mobilityRadius=0.
        gf = FirstPassageGreensFunction2D( self.getD() )
        self.domains = [ RadialDomain( self.getMobilityRadius(), gf ) ]


    def displacement( self, r ):
        x, y = randomVector2D( r )
        return x * self.surface.unitX + y * self.surface.unitY


    def __str__( self ):
        return 'PlanarSurfaceSingle' + Single.__str__( self )



'''
1 Particle inside a (cylindrical) shell on a CylindricalSurface. (Rods).

    * Particle coordinates on surface: z.
    * Domain: cartesian z.
    * Initial position: z=0.
    * Selected randomly when drawing displacement vector: none.
'''
class CylindricalSurfaceSingle( NonInteractionSingle ):
    def __init__( self, particle, reactionTypes ):
        NonInteractionSingle.__init__( self, particle, reactionTypes )

        # Heads up. The cylinder's *size*, not radius, is determined by 
        # getMinRadius(), because of redefinition of getRadius.
        self.shellList = [ Cylinder( particle.pos, self.surface.radius * SAFETY + 2 * particle.radius, self.surface.unitZ, self.getMinRadius() ) ]

        # Create a cartesian domain of size mobilityRadius=0.
        gf = FirstPassageGreensFunction1D( self.getD() )
        self.domains = [ CartesianDomain( 0, self.getMobilityRadius(), gf ) ]


    def displacement( self, z ):
        # z can be pos or min.
        return z * self.shellList[0].unitZ


    '''
    Overloaded methods getRadius and setRadius. Nice trick. property() needs 
    to be redefined as well.
    '''
    def getRadius( self ):
        # Heads up. Return cylinder's size.
        #raise RuntimeError, 'Todo'
        return self.shellList[0].size
    def setRadius( self, size ):
        #raise RuntimeError, 'Todo'
        assert size - self.getMinRadius() >= 0.0 # Still fine.
        # Heads up. A larger shell means a larger CylindricalSurfaceSingle's size.
        self.shellList[0].size = size
        self.domains[0].a = self.getMobilityRadius()
    radius = property( getRadius, setRadius )


    def __str__( self ):
        return 'CylindricalSurfaceSingle' + Single.__str__( self )


##############################################################################
'''
Interactions singles are used when a particle is close to a surface.

There are 2 types of InteractionSingles:

'''
class InteractionSingle( Single ):
    def __init__( self, particle, surface, reactionTypes, interactionType ):
        self.interactionSurface = surface
        self.interactionType = interactionType
        Single.__init__( self, particle, reactionTypes )


    # Interaction singles can not be updated, so no setters.
    def getPos( self ):
        return self.shellList[0].origin
    pos = property( getPos )


    def getRadius( self ):
        #raise RuntimeError, 'Todo'
        return self.shellList[0].radius
    radius = property( getRadius )


    def initialize( self, t ):
        self.lastTime = t


    def reset( self ):
        raise SystemError, 'Interaction singles should never be reset (reused)'


'''
1 Particle close to a PlanarSurface, inside a cylindrical shell placed on top of the surface.

    * Particle coordinates inside shell: r, theta, z.
    * Domains: radial r, cartesian z.
    * Initial position: r=0, z=z.
    * Selected randomly when drawing displacement vector: theta.
'''
class PlanarSurfaceInteraction( InteractionSingle ):
    def __init__( self, particle, surface, reactionTypes, interactionType, origin, radius, orientationVector, size, particleOffset, projectedPoint = None, sizeOfDomain = None ):
        InteractionSingle.__init__( self, particle, surface, reactionTypes, interactionType )

        self.shellList = [ Cylinder( origin, radius, orientationVector, size ) ]

        # Free diffusion in r direction.
        gfr = FirstPassageGreensFunction2D( self.getD() )
        rDomain = RadialDomain( self.getMobilityRadius(), gfr )

        '''
        Interaction possible in z direction.

        Heads up. sizeOfDomain is different from the size of the cylinder, 
        because the domain ends where the surface starts, while the 
        cylinder is continuing into the surface.
        '''
        gfz = FirstPassageGreensFunction1DRad( self.getD(), interactionType.k )
        # Cartesian domain of size L. Correction with getMinRadius().
        zDomain = CartesianDomain( particleOffset[1] - self.getMinRadius(), sizeOfDomain - 2*self.getMinRadius(), gfz )

        self.domains = [ rDomain, zDomain ]


    def drawNewPosition( self, dt ):
        r = self.domains[0].drawPosition( dt )
        # Cartesian domain returns displacement, not absolute position.
        z = self.domains[1].drawPosition( dt )
        x, y = randomVector2D( r )
        return self.pos + x * self.interactionSurface.unitX + y * self.interactionSurface.unitY + z * self.shellList[0].unitZ


    def __str__( self ):
        return 'PlanarSurfaceInteraction' + Single.__str__( self )


'''
1 Particle close to a CylindricalSurface, inside a cylindrical shell that surrounds the surface.

    * Particle coordinates inside shell: r, theta, z.
    * Domains: composite r-theta, cartesian z.
    * Initial position: r=r, theta=0, z=z.
    * Selected randomly when drawing displacement vector: none.
'''
class CylindricalSurfaceInteraction( InteractionSingle ):
    def __init__( self, particle, surface, reactionTypes, interactionType, origin, radius, orientationVector, size, particleOffset, projectedPoint, sizeOfDomain=None ):
        self.unitR = normalize( particle.pos - projectedPoint ) # Only needed for this type of Interaction.
        InteractionSingle.__init__( self, particle, surface, reactionTypes, interactionType )

        self.shellList = [ Cylinder( origin, radius, orientationVector, size ) ]

        # Interaction possible in r direction.
        self.pgf = FirstPassagePairGreensFunction2D( self.getD(), interactionType.k, surface.radius )
        rthetaDomain = CompositeDomain( surface.radius + self.getMinRadius(), particleOffset[0], radius - self.getMinRadius(), self.pgf )

        # Free diffusion in z direction.
        gfz = FirstPassageGreensFunction1D( self.getD() )
        # Cartesian domain of size L. Correction with getMinRadius().
        zDomain = CartesianDomain( particleOffset[1] - self.getMinRadius(), sizeOfDomain - 2*self.getMinRadius(), gfz )

        self.domains = [ rthetaDomain, zDomain ]


    def drawNewPosition( self, dt ):
        r, theta = self.domains[0].drawPosition( self.pgf, dt )
        # Cartesian domain returns displacement, not absolute position.
        z = self.domains[1].drawPosition( dt )
        # Calculate new position starting from origin.
        newVectorR = r * rotateVector( self.unitR, self.shellList[0].unitZ, theta )
        # Orientation matters (+ or -), so shellList[0].unitZ is used here 
        # instead of surface.unitZ.
        return self.pos + newVectorR + z * self.shellList[0].unitZ


    def __str__( self ):
        return 'CylindricalSurfaceInteraction' + Single.__str__( self )


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

