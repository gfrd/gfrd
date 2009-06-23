from _gfrd import EventType
from shape import Sphere
from utils import INF,NOWHERE
import math
import numpy

class Single( object ):
    def __init__( self, particle, shell, reactionTypes, domains ):
        self.particle = particle
        self.reactionTypes = reactionTypes
        self.domains = domains  # Can be composite.
        if reactionTypes:
            self.k_tot = sum( ( rt.k for rt in reactionTypes ) )
        else: 
            self.k_tot = 0
        self.shellList = [ shell ]
        self.eventType = None
        self.eventID = None
        self.multiplicity = 1


    def distanceTo( self, pos ):
        dists = numpy.zeros( len( self.shellList ) )
        for i, shell in enumerate( self.shellList ):
            dists[i] = shell.distanceTo( pos )
        return min( dists )


    def getD( self ):
        return self.particle.species.D


    def initialize( self, t ):
        self.reset()
        self.lastTime = t


    def burst( self, dt ):
        return self.toExternal( [ d.drawPosition( dt ) for d in self.domains ] )


    def propagate( self, dt ):
        return self.toExternal( [ d.drawPosition( dt ) for d in self.domains ] )


    # Determines a new event time and event type.
    def determineNextEvent( self, t ):
        self.dt, self.eventType = min(self.drawEscapeTime(), self.drawReactionTime()) 
        self.lastTime = t


    # Returns an (escapeTime, eventType)-tuple.
    # Handles also all interaction events.
    def drawEscapeTime( self ):
        if self.getD() == 0:
            return INF, EventType.ESCAPE
        else:
            t, d = min( (d.drawTime(), d) for d in self.domains )
            return t, d.drawEventType(t)


    # Returns a (reactionTime, eventType)-tuple.
    def drawReactionTime( self ):
        if self.k_tot == 0:
            dt = INF
        elif self.k_tot == INF:
            dt = 0.0
        else:
            dt = ( 1.0 / self.k_tot ) * math.log( 1.0 / numpy.random.uniform() )
        #return (dt, EventType.SINGLE_REACTION)
        # For now:
        return dt, EventType.REACTION


    # Copy pasted.
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
        return 'Shell' + str( self.particle )



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
