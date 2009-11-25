import numpy
import random
from _gfrd import EventType, FirstPassageGreensFunction1D
from log import *
from utils import *


'''
                radius of                                            radius of
            (surface) particle                                       particle
                |       |                                               |
                V       V                                               v
cartesian   | - - - |-------|-------------------|-------------------|-------|
domain            begin     0                   r0                  L      end
(also for          of                                                      of
 CoM and IV      shell                                                   shell
 1D pairs)




                                                                     radius of
                                                                     particle
                                                                        |
                                                                        V
radial                      |---------------------------------------|-------|
domain                    r0=0                                      a      end
                                                                           of
                                                                          shell



                    radii of particles                               radius of
                  (or particle+surface)                              particle
                            |                                           |
                            V                                           v
composite           |--------------|-------------|------------------|-------|
domain              0            sigma           r0                 a      end
(radial part)                                                              of
                                                                          shell
'''


class Domain( object ):
    """1 coordinate (2 in case of composite domain) of a position vector.

    The position of a particle inside a shell can be specified by a vector. 
    Which coordinates are used for this vector depends on the surface the 
    particle is on (if any). For example for a particle inside a sphere, not 
    on a surface, this is ( r, theta, phi ). But for a particle on a 
    CylindricalSurface, moving in 1 dimension, just ( z ) is sufficient.

    In the files single.py and pair.py it is defined which coordinates are 
    needed to define the position of the particle inside its shell for each 
    Single and Pair. Then, for each such coordinate, that can not be randomly 
    chosen, a domain is constructed.

    There are 3 types of domains:
        * CartesianDomain: for example the z-domain.
        * RadialDomain: for example the r-domain.
        * CompositeDomain: for example the r-theta-domain.

    All Greens Functions should be called from this file only.

    """
    def __init__( self, gf, a=None ):
        self.gf = gf         # Green's function.
        if a:
            self.a = a       # Outer radius.
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
    a = property( geta, seta )


class CartesianDomain( Domain ):
    """
    For example the z-domain. Extends from 0 to L.

    """
    def __init__( self, r0, L, gf ):
        Domain.__init__( self, gf )
        self.L = L
        self.r0 = r0         # Initial position


    def getr0( self ):
        return self._r0
    def setr0( self, r0 ):
        self._r0 = r0
        # Reset r0 of Greens' function here.
        self.gf.setr0( r0 )
    r0 = property( getr0, setr0 )


    def getL( self ):
        return self._L
    def setL( self, L ):
        self._L = L
        # Set L of Greens' function here.
        self.gf.setL( L )
    L = property( getL, setL )


    def drawTime( self ):
        assert self.escape == False
        try:
            rnd = numpy.random.uniform()
            log.debug( '        *Cartesian drawTime. ' ) #+ str( self.gf ) )
            dt = self.gf.drawTime( rnd )
        except Exception, e:
            raise RuntimeError( 'gf.drawTime() failed; %s; rnd = %g; r0 = %g; '
                                'L = %g' % ( str( e ), rnd, self.r0, self.L ) )
        return dt


    def drawEventType( self, dt ):
        """Return the type of event in this domain (which side) the particle 
        will leave.

        (!) Not a pure function.
        """

        # Set escape flag (can still be an interaction).
        self.escape = True
        try:
            log.debug( '        *Cartesian drawEventType. r0 = %.3g. '
                       'L = %.3g. dt = %.3g. ' %
                       ( self.r0, self.L, dt ) ) #+ str( self.gf ) )
            eventType = self.gf.drawEventType( numpy.random.uniform(), dt )
        except Exception, e:
            raise RuntimeError( 'gf.drawEventType() failed; %s; r0 = %g; '
                                'dt = %g; L = %g' %
                                ( str( e ), self.r0, dt, self.L ) )
        if eventType == EventType.REACTION:
            if( isinstance( self.gf, FirstPassageGreensFunction1D ) ):
                # If this is a 2xabsorbing greens function, eventType should 
                # always be ESCAPE.
                eventType = EventType.ESCAPE
            else:
                # Interaction.
                pass
            # Return displacement!
            self.newPos = 0 - self.r0
        elif eventType == EventType.ESCAPE:
            # Return displacement!
            self.newPos = self.L - self.r0
        return eventType

    
    def drawPosition( self, dt ):
        """Returns displacement, not absolute position in this domain!
        (!) Not a pure function.

        """
        if self.escape:
            # Escape through this domain. We already know which side.
            self.escape = False
            return self.newPos

        try:
            rnd = numpy.random.uniform()
            log.debug( '        *Cartesian drawR. ' ) #+ str( self.gf ) )
            r = self.gf.drawR( rnd, dt  )
        except Exception, e:
            raise RuntimeError( 'gf.drawR failed; %s; rnd = %g, dt = %g,'
                                'r0 = %g, L = %g' %
                                ( str( e ), rnd, dt, self.r0, self.L ) )

        # Return displacement.
        return r - self.r0


class RadialDomain( Domain ):
    """For example the r-domain. 

    The initial position of the particle is always at r = 0, so theta can and 
    should be choosen at random.

    """
    def __init__( self, a, gf ):
        Domain.__init__( self, gf, a )


    def drawTime( self ):
        try:
            rnd = numpy.random.uniform()
            log.debug( '        *Radial drawTime. ' ) #+ str( self.gf ) )
            dt = self.gf.drawTime( rnd )
        except Exception, e:
            raise RuntimeError( 'gf.drawTime() failed; %s; rnd = %g; a = %g' %
                                ( str( e ), rnd, self.a ) )
        return dt


    def drawEventType( self, _ ):
        """(!) Not a pure function.

        """
        self.escape = True
        # Always return escape.
        return EventType.ESCAPE


    def drawPosition( self, dt ):
        """(!) Not a pure function.

        """
        if self.escape:
            # Escape through this domain. We already know the new r.
            self.escape = False
            return self.a
        try:
            rnd = numpy.random.uniform()
            log.debug( '        *Radial drawR. ' ) #+ str( self.gf ) )
            r = self.gf.drawR( rnd, dt  )
        except Exception, e:
            raise RuntimeError( 'gf.drawR failed; %s; rnd = %g, dt = %g, '
                                'a = %g' % ( str( e ), rnd, dt, self.a ) )

        return r


class CompositeDomain( Domain ):
    """For example the ( r, theta ) domain used with PairGreensFunctions (3D 
    as well as 2D).

    """
    def __init__( self, sigma, r0, a, gf ):
        self.sigma = sigma              # Inner radius.
        self.r0 = r0                    # Starting position.
        Domain.__init__( self, gf, a )  # Greens function holds D, sigma, k.


    def drawTime( self ):
        assert self.escape == False
        try:
            rnd = numpy.random.uniform()
            log.debug( '        *Radial2D drawTime. ' ) #+ str( self.gf ) )
            dt = self.gf.drawTime( rnd, self.r0 )
        except Exception, e:
            raise RuntimeError( 'gf.drawTime() failed; %s; rnd = %g, '
                                'sigma = %g, r0 = %g, a = %g' %
                                ( str( e ), rnd, self.sigma, self.r0, self.a ) )
        self.check_dt = dt
        self.check_gf = self.gf
        return dt


    def drawEventType( self, dt ):
        """(!) Not a pure function.

        """
        try:
            rnd = numpy.random.uniform()
            log.debug( '        *Radial2D drawEventType. ' ) #+ str( self.gf ) )
            eventType = self.gf.drawEventType( rnd, self.r0, dt )
        except Exception, e:
            raise RuntimeError( 'gf.drawEventType() failed; %s; sigma = %g;'
                                'r0 = %g; a = %g; dt = %g' %
                                ( str( e ), self.sigma, self.r0, self.a, dt ) )
        if eventType == EventType.ESCAPE:
            self.escape = True
        return eventType     # 0 (REACTION) or 1 (ESCAPE r)


    def drawPosition( self, gf, dt ):
        """(!) Not a pure function.

        """
        if self.escape:
            self.escape = False
            # Escape through this domain. We already know the new r.
            r = self.a
        else:
            r = self.drawR_pair( gf, dt )
        theta = self.drawTheta_pair( gf, r, dt )
        return r, theta


    def drawR_pair( self, gf, dt ):
        """Draw r for the pair inter-particle vector.

        """
        try:
            rnd = numpy.random.uniform()
            log.debug( '        *Radial2D drawR_pair. ' ) #+ str( self.gf ) )
            r = gf.drawR( rnd, self.r0, dt )
            # redraw; shouldn't happen often
            while r >= self.a or r <= self.sigma: 
                #log.info( '    drawR_pair: redraw' )
                #self.sim.rejectedMoves += 1  #FIXME:
                rnd = numpy.random.uniform()
                r = gf.drawR( rnd, self.r0, dt )
        except Exception, e:
            raise RuntimeError( 'gf.drawR_pair() failed; %s; rnd = %g, '
                                'sigma = %g, r0 = %g, a = %g, dt = %g' % 
                                ( str( e ), rnd, self.sigma, self.r0, self.a, 
                                  dt ) )
        return r


    def drawTheta_pair( self, gf, r, dt ):
        """Draw theta for the pair inter-particle vector.

        """
        try:
            rnd = numpy.random.uniform()
            log.debug( '        *Radial2D drawTheta_pair. ' )#+ str( self.gf ) )
            theta = gf.drawTheta( rnd, r, self.r0, dt )
        except Exception, e:
            raise RuntimeError( 'gf.drawTheta() failed; %s; rnd = %g, r = %g, '
                                'sigma = %g, r0 = %g, a = %g, dt = %g' %
                                ( str( e ), rnd, r, self.sigma, self.r0, 
                                  self.a, dt  ) )

        # Heads up. For cylinders theta should be between [ -pi, pi ]. For 
        # spheres it doesn't matter.
        return random.choice( [ -1, 1 ] ) * theta


