import numpy
import random
from _gfrd import EventType


class Domain( object ):
    def __init__( self, a, gf ):
        self.gf = gf         # Green's function.
        self.a = a           # Out er radius.
        self.escape = False  # Flag for if there is an escape in this domain 
                             # later. So we don't have to draw a position for 
                             # this domain in case of a propagate. 
    

    def drawTime( self ):
        # Todo. This only works this simple if all parameters have been set in 
        # greens function, like gf.seta(), but also gf.ka() and gf.kb().
	try:
            rnd = numpy.random.uniform()
	    dt = self.gf.drawTime( rnd )
	except Exception, e:
	    raise Exception, 'gf.drawTime() failed; %s; rnd=%g, %s' %\
		( str( e ), rnd, self.gf.dump() )
	return dt


    def geta( self ):
        return self._a
    def seta( self, a ):
        self._a = a
        self.gf.seta( a )
    a = property(geta, seta)


    def anyGreensFunctionDrawR( self, dt ):
        try:
            rnd = numpy.random.uniform()
            r = self.gf.drawR( rnd , dt )
        except Exception, e:
            raise Exception, 'gf.drawR failed; %s; rnd=%g, t=%g, %s' %\
                ( str( e ), rnd, dt, self.gf.dump() )
        return r



'''
For example the x-domain.
'''
class CartesianDomain1D( Domain ):
    def __init__( self, r0, (kInner, kOuter), a, gf ):
        self.r0 = r0         # Initial position
        self.kInner = kInner # Inner interaction rate.
        self.kOuter = kOuter # Outer interaction rate.
        self.newPos = a      # Position in this domain where particle will end 
                             # up after an escape.
        Domain.__init__( self, a, gf )


    '''
    Returns the type of event in this domain (which side) the particle will 
    leave.
    TODO
    (!) Not a pure function.
    '''
    def drawEventType( self, dt ):
        self.escape = True
        # Todo after here. Flux stuff.
        choice = random.choice( [-1,1] )
        self.newPos = self.a * choice
        if choice == -1:
            return EventType.REACTION
        else:
            return EventType.ESCAPE

        try:
            eventType = self.gf.drawEventType( numpy.random.uniform(), self.r0, 
                self.kInner, self.kOuter, dt )
	except Exception, e:
            raise Exception, 'gf.drawEventType() failed; %s; r0=%g, %s' %\
                ( str( e ), self.r0, self.pgf.dump() )
        if eventType == EventType.REACTION:
            # Interaction.
            self.newPos = - self.a
        elif eventType == EventType.ESCAPE:
            # Escape.
            self.newPos = self.a
        return eventType

    
    # (!) Not a pure function.
    def drawPosition( self, dt ):
        if self.escape:
            # Escape through this domain. We already know which side.
            self.escape = False
            return self.newPos
        r = self.anyGreensFunctionDrawR( dt )
        return r * random.choice( [-1,1] ) # Cartesian domain. Move either left or right.


class RadialDomain1D( Domain ):
    def __init__( self, a, gf ):
        Domain.__init__( self, a, gf )


    # (!) Not a pure function.
    def drawEventType( self, _ ):
	self.escape = True
	return EventType.ESCAPE


    # (!) Not a pure function.
    def drawPosition( self, dt ):
        if self.escape:
            # Escape through this domain. We already know the new r.
            self.escape = False
            return self.a
	return self.anyGreensFunctionDrawR( dt )


'''
For example the (r-theta) domain.
'''
class RadialDomain2D( Domain ):
    def __init__( self, sigma, r0, a, gf ):
        self.sigma = sigma              # Inner radius.
        self.r0 = r0                    # Starting position.
        Domain.__init__( self, a, gf )  # Greens function holds D, sigma, k.

    # Overloaded method.
    def drawTime( self ):
        assert self.escape == False
        # Todo. This only works this simple if all parameters have been set in 
        # greens function, like gf.seta(), but also gf.ka() and gf.kb().
        print 'r0 = ', self.r0
	try:
            rnd = numpy.random.uniform()
	    dt = self.gf.drawTime( rnd, self.r0 )
	except Exception, e:
	    raise Exception, 'gf.drawTime() failed; %s; rnd=%g, %s' %\
		( str( e ), rnd, self.gf.dump() )
        self.check_dt = dt
        self.check_gf = self.gf
	return dt


    # (!) Not a pure function.
    def drawEventType( self, dt ):
        try:
            rnd = numpy.random.uniform()
            eventType = self.gf.drawEventType( rnd, self.r0, dt )
        except Exception, e:
            raise Exception,\
                'gf.drawEventType() failed; %s; r0=%g, %s' %\
                ( str( e ), self.r0, self.gf.dump() )
        if eventType == EventType.ESCAPE:
            self.escape = True
        return eventType     # 0 (REACTION) or 1 (ESCAPE r)


    # (!) Not a pure function.
    def drawPosition( self, gf, dt ):
        if self.escape:
            self.escape = False
            # Escape through this domain. We already know the new r.
            r = self.a
        else:
            r = self.drawR_pair( gf, dt )
        theta = self.drawTheta_pair( gf, r, dt )
        return (r, theta)


    # Draw r for the pair inter-particle vector.
    def drawR_pair( self, gf, dt ):
        try:
            rnd = numpy.random.uniform()
            r = gf.drawR( rnd, self.r0, dt )
            # redraw; shouldn't happen often
            while r >= self.a or r <= self.sigma: 
                #log.info( 'drawR_pair: redraw' )
                #self.sim.rejectedMoves += 1  #FIXME:
                rnd = numpy.random.uniform()
                r = gf.drawR( rnd, self.r0, dt )
        except Exception, e:
            raise Exception,\
                'gf.drawR_pair() failed; %s; rnd= %g, r0= %g, dt= %g, %s' %\
                ( str( e ), rnd, self.r0, dt, gf.dump() )
        return r


    # Draw theta for the pair inter-particle vector.
    def drawTheta_pair( self, gf, r, dt ):
        try:
            rnd = numpy.random.uniform()
            theta = gf.drawTheta( rnd, r, self.r0, dt )
        except Exception, e:
            raise Exception,\
                'gf.drawTheta() failed; %s; rnd= %g, r= %g, r0= %g, dt=%g, %s' %\
                ( str( e ), rnd, r, self.r0, dt, gf.dump() )

        '''
        Heads up. For cylinders theta should be between [-pi,pi]. For spheres 
        it doesn't matter.
        '''
        return random.choice( [-1,1] ) * theta


