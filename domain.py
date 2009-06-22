import numpy
from _gfrd import EventType

class Domain( object ):
    # Todo: bundle more stuff from 3 classes below.
    def __init__( self ):
        pass


# For example the x-domain, or the r-domain (x0=0).
class SimpleDomain( Domain ):
    def __init__( self, x0, a, (k0, ka), gf ):
        self.x0 = x0        # Initial position.
        self.a = a          # Length or radius, a.k.a. size.
                            # Note: this should be the mobilityRadius.
        self.k0 = k0        # Inner interaction rate.
                            # Todo: take into account mobility radius here as 
                            # well. This is important. So subtract particle 
                            # radius twice or something.
        self.ka = ka        # Outer interaction rate.
        self.gf = gf        # Green's function.
        self.active = True  # Flag for if there is an escape through this domain.
        self.newPos = a     # Todo. FixMe.

    def drawPosition( self, dt ):
        # Escape through this domain. We already know the new value.
        if self.active:
            self.active = False
            return self.newPos

        rnd = numpy.random.uniform()
        self.gf.seta( self.a )
        if self.x0 == 0 and self.k0 == None and self.ka == 0:
            # FreeSingle.
            try:
                r = self.gf.drawR( rnd , dt )
            except Exception, e:
                raise Exception, 'gf.drawR failed; %s; rnd=%g, t=%g, %s' %\
                    ( str( e ), rnd, dt, self.gf.dump() )
            return r
        else:
            # General case.
            return self.gf.drawR( numpy.random.uniform(), self.x0, self.k0, self.ka, dt )

    def drawTime( self ):
        self.gf.seta( self.a )
        rnd = numpy.random.uniform()
        if self.x0 == 0 and self.k0 == None and self.ka == 0:
            # FreeSingle.
            try:
                dt = self.gf.drawTime( rnd )
            except Exception, e:
                raise Exception, 'gf.drawTime() failed; %s; rnd=%g, %s' %\
                    ( str( e ), rnd, self.gf.dump() )
            return dt
        else:
            # General case.
            return self.gf.drawTime( numpy.random.uniform(), self.x0, self.k0, self.ka )

    # Returns the type of event and the position in this domain (which side)    
    # the particle will leave. The first return value should be either:
    # * EventType.INTERACTION (in case radiating boundary condition at r=0) = 
    # REACTION (in case pair) = ESCAPE_0 (in case absorbing boundary condition 
    # at r=0)
    # * or EventType.ESCAPE
    def drawEventType( self, t ):
        self.gf.seta( self.a )
        if self.x0 == 0 and self.k0 == None and self.ka == 0:
            # FreeSingle.
            self.active = True
            self.newPos = self.a
            return EventType.ESCAPE
        else:
            # General case:
            # Todo, set:
            # self.newPos =
            # self.active =
            return self.gf.drawEventType( numpy.random.uniform(), self.x0, self.k0, self.ka, t )



# For example the (r-theta) domain.
class CompositeDomain( Domain ):
    def __init__( self, r0, sigma, a, (ksigma, ka), gf ):
        self.r0 = r0            # Initial position.
        self.sigma = sigma      # Inner radius.
        self.a = a              # Outer radius.
        self.ksigma = ksigma    # Inner interaction rate.
        self.ka = ka            # Outer interaction rate.
        self.gf = gf            # Green's function.

    def drawPosition( self, dt ):
	# Todo: exception handling.
        gf = self.choosePairGreensFunction( self.r0, t )
        gf.seta( self.a )
        r = self.gf.drawR( numpy.random.uniform(), self.r0, dt )
        theta = self.gf.drawTheta( numpy.random.uniform(), r, self.r0, dt )
        return (r, theta)

    def drawTime( self ):
        self.seta( self.a )
	# Todo: exception handling.
        return self.gf.drawTime( numpy.random.uniform(), self.r0 )







################## GRAVEYARD #####################


# For example the r-domain.
class RadialDomain( Domain ):
    def __init__( self, a, gf ): 
        self.a = a
        self.gf = gf

    def drawPosition( dt ):
        self.gf.seta( self.a )
        return self.gf.drawX( dt )


class Boundary( object ):
    def __init__( self, r, k=None, surface=None ):
        self.r = r          # For example: r=a.    
        self.k = k          # Reaction rate.
        self.s = surface    # Target surface.

