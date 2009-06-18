import numpy

from utils import *
from shell import *
from _gfrd import *
from gfrdbase import log


class Single( object ):
    global log

    def __init__( self, particle, reactiontypes ):

        self.multiplicity = 1

        self.particle = particle
        self.reactiontypes = reactiontypes

        self.k_tot = 0

        self.lastTime = 0.0
        ### self.dt will be firstPassageTime or reactionTime when
        ### updateNextEvent() is called.
        self.dt = 0.0
        self.eventType = None

        ### So shell radius is particle radius.
        ### This is enlarged once updateSingle() is called.
        ### Note it's possible for a Single/Pair/Multi to have multiple
        ### shells.
        self.shellList = [ Shell( self.particle.pos, self.getMinRadius() ), ]

        self.eventID = None

        ### Green's function for particle inside absorbing sphere. C++.
        self.gf = FirstPassageGreensFunction( particle.species.D )

        ### Total of all reaction constants.
        self.updatek_tot()


    def getD( self ):

        return self.particle.species.D


    def getPos( self ):

        return self.shellList[0].pos


    def setPos( self, pos ):
        self.shellList[0].pos = pos
        self.particle.pos = pos

    pos = property( getPos, setPos )


    def setRadius( self, radius ):

        assert radius - self.getMinRadius() >= 0.0

        self.shellList[0].radius = radius


    '''
    A radius of a Single is the distance from the current position
    of the particle to the shell of the Single.   The shell defines the
    farthest point in space that it can occupy when it makes the maximum
    displacement.
    '''

    def getRadius( self ):

        return self.shellList[0].radius

    radius = property( getRadius, setRadius )


    ### This will also be the radius of a shell after initialization.
    def getMinRadius( self ):

        return self.particle.species.radius



    '''
    Initialize this Single.

    The radius (shell size) is shrunken to the radius of the particle
    it represents.   
    self.lastTime is reset to the current time, and self.dt
    is set to zero.
    '''

    def initialize( self, t ):

        self.reset()
        self.lastTime = t



    '''
    A mobility radius indicates the maximum displacement this single
    particle can make.

    Mobility radius of a particle is calculated as follows;

    mobility radius = Single radius - Particle radius.

    '''
    
    def getMobilityRadius( self ):

        #return self.radius - ( self.getMinRadius() * 2 )
        return self.radius - self.getMinRadius()




    '''
    Reset the Single.

    Radius (shell size) is shrunken to the actual radius of the particle.
    self.dt is reset to 0.0.  Do not forget to reschedule this Single
    after calling this method.
    '''

    def reset( self ):

        self.setRadius( self.getMinRadius() )
        self.dt = 0.0
        self.eventType = EventType.ESCAPE


    def isReset( self ):

        return self.radius == self.getMinRadius() and self.dt == 0.0\
               and self.eventType == EventType.ESCAPE
        

    ### This is not free diffusion, but diffusion with absorbing boundary
    ### conditions.
    def drawR( self, dt ):
        assert dt >= 0

        rnd = numpy.random.uniform()
        ### Not needed, already in self.gf.drawR()?
        self.gf.seta( self.getMobilityRadius() )

        try:
            r = self.gf.drawR( rnd , dt )
        except Exception, e:
            raise Exception, 'gf.drawR failed; %s; rnd=%g, t=%g, %s' %\
                ( str( e ), rnd, dt, self.gf.dump() )

        return r



    ### Called by EGFRDSimulator::updateSingle() and 
    ### EGFRDSimulator::fireSingle() (for immobile case only).
    def determineNextEvent( self, t ):

        if self.getD() == 0:
            firstPassageTime = INF
        else:
            firstPassageTime = self.drawEscapeTime()
            
        reactionTime = self.drawReactionTime()

        if firstPassageTime <= reactionTime:
            ### ! Here Single::dt is set.
            self.dt = firstPassageTime
            self.eventType = EventType.ESCAPE
        else:
            self.dt = reactionTime
            self.eventType = EventType.REACTION

        self.lastTime = t


    ### Called by Single::determineNextEvent().
    def drawReactionTime( self ):
        
        if self.k_tot == 0:
            return INF

        if self.k_tot == INF:
            return 0.0

        rnd = numpy.random.uniform()
        ### Notes December 1 2008.
        dt = ( 1.0 / self.k_tot ) * math.log( 1.0 / rnd )

        return dt


    ### Called by Single::determineNextEvent().
    def drawEscapeTime( self ):
        
        rnd = numpy.random.uniform()

        self.gf.seta( self.getMobilityRadius() )

        try:
            dt = self.gf.drawTime( rnd )
        except Exception, e:
            raise Exception, 'gf.drawTime() failed; %s; rnd=%g, %s' %\
                ( str( e ), rnd, self.gf.dump() )

        return dt


    def updatek_tot( self ):

        self.k_tot = 0

        if not self.reactiontypes:
            return

        for rt in self.reactiontypes:
            self.k_tot += rt.k


    ### Called by fireSingle().
    def drawReactionType( self ):

        k_array = [ rt.k for rt in self.reactiontypes ]
        k_array = numpy.add.accumulate( k_array )
        ### This should be equal to self.k_tot? No, it's a cumulative addition
        ### me thinks.
        k_max = k_array[-1]

        rnd = numpy.random.uniform()
        i = numpy.searchsorted( k_array, rnd * k_max )

        return self.reactiontypes[i]


    def check( self ):
        pass

    def __str__( self ):
        return 'Single' + str( self.particle )




class DummySingle( object ):
    def __init__( self ):
        self.multiplicity = 1

        self.radius = 0.0
        self.shellList = [ Shell( NOWHERE, 0.0 ), ]


    def getMinRadius( self ):
        return 0.0

    def getD( self ):
        return 0.0

    def getPos( self ):
        return NOWHERE

    pos = property( getPos )


    def __str__( self ):
        return 'DummySingle()'



