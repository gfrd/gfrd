import numpy
import random
from _gfrd import EventType
#from gfrdbase import log


class Domain( object ):
    def __init__( self, a, gf ):
        self.gf = gf         # Green's function.
        self.a = a           # Out er radius.
        self.escape = False  # Flag for if there is an escape in this domain 
                             # later. So we don't have to draw a position for 
                             # this domain in case of a propagate. By default 
                             # all escape flags have to be set to False, so 
                             # drawPosition works properly.
    

    def geta( self ):
        return self._a
    def seta( self, a ):
        self._a = a
        # Set a of Greens' function here.
        self.gf.seta( a )
    a = property(geta, seta)


'''
For example the z-domain.
'''
class CartesianDomain1D( Domain ):
    def __init__( self, r0, a, gf ):
        self.r0 = r0         # Initial position
        Domain.__init__( self, a, gf )


    def drawTime( self ):
        assert self.escape == False
	try:
            rnd = numpy.random.uniform()
            #print( '\tDebug. Cartesian drawTime.' )
	    dt = self.gf.drawTime( rnd, self.r0 )
	except Exception, e:
	    raise Exception, 'gf.drawTime() failed; %s; rnd=%g; r0=%g; a=%g' %\
		( str( e ), rnd, self.r0, self.a )
	return dt


    '''
    Returns the type of event in this domain (which side) the particle will 
    leave.
    (!) Not a pure function.
    '''
    def drawEventType( self, dt ):
        # Set escape flag.
        self.escape = True
        try:
            #print( '\tDebug. Cartesian drawEventType.' )
            eventType = self.gf.drawEventType( numpy.random.uniform(), self.r0, dt )
	except Exception, e:
            raise Exception, 'gf.drawEventType() failed; %s; r0=%g; dt=%g; a=%g' %\
                ( str( e ), self.r0, self.dt, self.a )
        if eventType == EventType.REACTION:
            # Interaction.
            self.newPos = - self.a
        elif eventType == EventType.ESCAPE:
            # Escape.
            self.newPos = self.a
        return eventType

    
    # (!) Not a pure function.
    # Returns displacement, not absolute position in this domain!
    def drawPosition( self, dt ):
        if self.escape:
            # Escape through this domain. We already know which side.
            self.escape = False
            return self.newPos

        try:
            rnd = numpy.random.uniform()
            #print( '\tDebug. Cartesian drawR.' )
            r = self.gf.drawR( rnd, self.r0, dt  )
        except Exception, e:
            raise Exception, 'gf.drawR failed; %s; rnd=%g, dt=%g, r0=%g, a=%g' %\
                ( str( e ), rnd, dt, self.r0, self.a )

        # Return displacement.
        return r - self.r0


class RadialDomain1D( Domain ):
    def __init__( self, a, gf ):
        Domain.__init__( self, a, gf )


    def drawTime( self ):
	try:
            rnd = numpy.random.uniform()
            #print( '\tDebug. Radial drawTime.' )
	    dt = self.gf.drawTime( rnd )
	except Exception, e:
	    raise Exception, 'gf.drawTime() failed; %s; rnd=%g; a=%g' %\
		( str( e ), rnd, self.a )
	return dt


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
        try:
            rnd = numpy.random.uniform()
            #print( '\tDebug. Radial drawR.' )
            r = self.gf.drawR( rnd, dt  )
        except Exception, e:
            raise Exception, 'gf.drawR failed; %s; rnd=%g, dt=%g, a=%g' %\
                ( str( e ), rnd, dt, self.a )

	return r


'''
For example the (r-theta) domain.
'''
class RadialDomain2D( Domain ):
    def __init__( self, sigma, r0, a, gf ):
        self.sigma = sigma              # Inner radius.
        self.r0 = r0                    # Starting position.
        Domain.__init__( self, a, gf )  # Greens function holds D, sigma, k.


    def drawTime( self ):
        assert self.escape == False
	try:
            rnd = numpy.random.uniform()
            #print( '\tDebug. Radial2D drawTime.' )
	    dt = self.gf.drawTime( rnd, self.r0 )
	except Exception, e:
	    raise Exception, 'gf.drawTime() failed; %s; rnd=%g, sigma=%g, r0=%g, a=%g' %\
		( str( e ), rnd, self.sigma, self.r0, self.a )
        self.check_dt = dt
        self.check_gf = self.gf
	return dt


    # (!) Not a pure function.
    def drawEventType( self, dt ):
        try:
            rnd = numpy.random.uniform()
            #print( '\tDebug. Radial2D drawEventType.' )
            eventType = self.gf.drawEventType( rnd, self.r0, dt )
        except Exception, e:
            raise Exception,\
                'gf.drawEventType() failed; %s; sigma=%g; r0=%g; a=%g; dt=%g' %\
                ( str( e ), self.sigma, self.r0, self.a, self.dt )
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
            #print( '\tDebug. Radial2D drawR_pair.' )
            r = gf.drawR( rnd, self.r0, dt )
            # redraw; shouldn't happen often
            while r >= self.a or r <= self.sigma: 
                #log.info( 'drawR_pair: redraw' )
                #self.sim.rejectedMoves += 1  #FIXME:
                rnd = numpy.random.uniform()
                r = gf.drawR( rnd, self.r0, dt )
        except Exception, e:
            raise Exception,\
                'gf.drawR_pair() failed; %s; rnd= %g, sigma=%g, r0= %g, a=%g, dt= %g' %\
                ( str( e ), rnd, self.sigma, self.r0, self.a, dt )
        return r


    # Draw theta for the pair inter-particle vector.
    def drawTheta_pair( self, gf, r, dt ):
        try:
            rnd = numpy.random.uniform()
        #print( '\tDebug. Radial2D drawTheta_pair.' )
            theta = gf.drawTheta( rnd, r, self.r0, dt )
        except Exception, e:
            raise Exception,\
                'gf.drawTheta() failed; %s; rnd= %g, r=%g, sigma=%g, r0=%g, a=%g, dt=%g' %\
                ( str( e ), rnd, r, self.sigma, self.r0, self.a, dt  )

        '''
        Heads up. For cylinders theta should be between [-pi,pi]. For spheres 
        it doesn't matter.
        '''
        return random.choice( [-1,1] ) * theta


