from utils import *
#from utils import cyclicTranspose, math, numpy, vectorAngleAgainstZAxis, crossproductAgainstZAxis, normalize, rotateVector, distance_Simple
from shell import *
from _gfrd import *
from gfrdbase import log

SAFETY = 1.0 + 1e-5


'''
Just a free func ver of Pair.getCoM().
'''

### Called by formPair().
def calculatePairCoM( pos1, pos2, D1, D2, worldSize ):

    pos2t = cyclicTranspose( pos2, pos1, worldSize )

    return ( ( D2 * pos1 + D1 * pos2t ) / ( D1 + D2 ) ) % worldSize


class Pair( object ):
    global log
    
    ### How does this work exactly?
    # CUTOFF_FACTOR is a threshold to choose between the real and approximate
    # Green's functions.
    # H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    # 5.6: ~1e-8, 6.0: ~1e-9
    CUTOFF_FACTOR = 5.6

    # The distance function is passed to be able to take the boundary 
    # conditions into account.
    def __init__( self, single1, single2, rt, distFunc, worldSize ):

        self.multiplicity = 2

        ### So a pair is a combination of 2 singles! (individual shells are
        ### preserved I suppose).
        # Order single1 and single2 so that D1 < D2.
        if single1.particle.species.D <= single2.particle.species.D:
            self.single1, self.single2 = single1, single2 
        else:
            self.single1, self.single2 = single2, single1 

        self.rt = rt

        self.distance = distFunc
        self.worldSize = worldSize
        
        particle1 = self.single1.particle
        particle2 = self.single2.particle

        self.D1, self.D2 = particle1.species.D, particle2.species.D

        self.D_tot = self.D1 + self.D2
        self.D_geom = math.sqrt( self.D1 * self.D2 )  # geometric mean

        #self.minRadius = max( particle1.species.radius,
        #                      particle2.species.radius )
        self.sigma = particle1.species.radius + particle2.species.radius

        ### Green's function for particle inside absorbing sphere. C++.
        self.sgf = FirstPassageGreensFunction( self.D_geom )
        ### Green's function for pair inside absorbing sphere? C++.
        self.pgf = FirstPassagePairGreensFunction( self.D_tot, 
                                                   rt.k, self.sigma )

        self.eventID = None

        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None

        self.shellList = [ Shell( self.getCoM(), self.getMinRadius() ), ]



    def __del__( self ):
        log.debug( 'del %s' % str( self ) )

    def initialize( self, t ):

        self.lastTime = t
        self.radius = self.getMinRadius()
        self.dt = 0
        self.eventType = None

    ### So the first shell determines the position of the pair.
    def setPos( self, pos ):
        self.shellList[0].origin = pos

    def getPos( self ):

        return self.shellList[0].origin

    pos = property( getPos, setPos )

    def getD( self ):

        return self.D_tot #FIXME: is this correct?

    def setRadius( self, radius ):

        assert radius - self.getMinRadius() >= 0.0
        self.shellList[0].size = radius


    def getRadius( self ):

        return self.shellList[0].size

    radius = property( getRadius, setRadius )


    '''
    This method returns the radius from its CoM that this Pair must reserve
    to remain mobile.
    '''

    ### Isn't this unneccesary, setting a minimum radius at instantiation?
    ### Will be overwritten anyway by a call to setRadius().
    def getMinRadius( self ):

        pairDistance = self.distance( self.single1.pos,
                                      self.single2.pos )
        minRadius = max( pairDistance * self.D1 /
                         self.D_tot + self.single1.getMinSize(),
                         pairDistance * self.D2 /
                         self.D_tot + self.single2.getMinSize() )
        return minRadius


    '''
    Calculate and return the "Center of Mass" (== CoM) of this pair.
    '''

    ### Used in constructor.
    def getCoM( self ):

        particle1 = self.single1.particle
        particle2 = self.single2.particle
        
        pos1 = particle1.pos
        pos2 = particle2.pos

        ### ?
        pos2t = cyclicTranspose( pos2, pos1, self.worldSize ) #FIXME:
        
        com = ( pos1 * self.D2 + pos2t * self.D1 ) / self.D_tot

        return com % self.worldSize


    ### Called by drawR_pair(), and drawTheta_pair().
    ### Decides if an approximate green's function is used
    ### ro is interparticle vector.
    def choosePairGreensFunction( self, r0, t ):

        ### sigma is particle1.radius+particle2.radius.
        distanceFromSigma = r0 - self.sigma
        distanceFromShell = self.a_r - r0

        thresholdDistance = Pair.CUTOFF_FACTOR * \
            math.sqrt( 6.0 * self.D_tot * t )

        if distanceFromSigma < thresholdDistance:
        
            ### Is it fair to compare distanceFromSigma as well as
            ### distanceFromShell to the same threshold value?
            ### I guess it is. How close are shells to each other.
            if distanceFromShell < thresholdDistance:
                # near both a and sigma;
                # use FirstPassagePairGreensFunction
                log.debug( 'GF: normal' )
                return self.pgf
            else:
                # near sigma; use BasicPairGreensFunction
                log.debug( 'GF: only sigma' )
                pgf = BasicPairGreensFunction( self.D_tot, self.rt.k, 
                                               self.sigma )
                return pgf
                #return self.pgf
        else:
            if distanceFromShell < thresholdDistance:
                # near a;
                log.debug( 'GF: only a' )
                pgf = FirstPassageNoCollisionPairGreensFunction( self.D_tot )
                return pgf
                
            else:
                # distant from both a and sigma; 
                log.debug( 'GF: free' )
                pgf = FreePairGreensFunction( self.D_tot )
                return pgf


    '''
    Calculate new positions of the pair particles using
    a new center-of-mass, a new inter-particle vector, and
    an old inter-particle vector.

    '''

    ### Called by firePair() and breakUpPair().
    ### newInterParticle is still something like [r, dtheta, dphi]. So here we
    ### find out how much we have to rotate it (taking into account the
    ### oldInterParticle) to get it's final orientation.
    def newPositions( self, CoM, newInterParticle, oldInterParticle ):

        #FIXME: need better handling of angles near zero and pi.

        # I rotate the new interparticle vector along the
        # rotation axis that is perpendicular to both the
        # z-axis and the original interparticle vector for
        # the angle between these.
        
        # the rotation axis is a normalized cross product of
        # the z-axis and the original vector.
        # rotationAxis = crossproduct( [ 0,0,1 ], interParticle )

        angle = vectorAngleAgainstZAxis( oldInterParticle )
        if angle % numpy.pi != 0.0:
            rotationAxis = crossproductAgainstZAxis( oldInterParticle )
            rotationAxis = normalize( rotationAxis )
            rotated = rotateVector( newInterParticle,
                                    rotationAxis,
                                    angle )
        elif angle == 0.0:
            rotated = newInterParticle
        else:
            rotated = numpy.array( [ newInterParticle[0], newInterParticle[1],
                                     - newInterParticle[2] ] )
            #rotated = newInterParticle * numpy.array( [ 1, 1, -1 ] )

        newpos1 = CoM - rotated * ( self.D1 / self.D_tot )
        newpos2 = CoM + rotated * ( self.D2 / self.D_tot )

        return newpos1, newpos2
        

    ### Called by formPair().
    def determineNextEvent( self, t ):
        
        self.lastTime = t

        single1 = self.single1
        single2 = self.single2
        radius1 = single1.particle.species.radius
        radius2 = single2.particle.species.radius

        D1 = self.D1
        D2 = self.D2

        D1_factor = D1 / self.D_tot
        D2_factor = D2 / self.D_tot

        shellSize = self.radius / SAFETY  # FIXME:

        sqrtD_tot = math.sqrt( self.D_tot )
        sqrtD_geom = math.sqrt( self.D_geom )

        r0 = self.distance( single1.pos, single2.pos )

        assert r0 >= self.sigma, \
            '%s;  r0 %g < sigma %g' % ( self, r0, self.sigma )

        # equalize expected mean t_r and t_R.

        ### Determine a_r and a_R?
        qrrtD1D25 = ( D1    * D2**5 ) ** 0.25
        qrrtD15D2 = ( D1**5 * D2 ) ** 0.25

        if qrrtD15D2 * r0 + ( qrrtD15D2 + qrrtD1D25 ) * radius1 \
                + D1 * ( sqrtD_tot * ( shellSize - radius2 ) 
                         - sqrtD_geom * radius2 )\
                - D2 * ( sqrtD_geom * r0 + sqrtD_tot * 
                         ( shellSize - radius1 ) )\
                         - qrrtD1D25 * radius2 >= 0:

            den1 = qrrtD1D25 + D1 * ( sqrtD_geom + sqrtD_tot )

            a_R_1 = sqrtD_geom * ( D2 * ( shellSize - radius1) + 
                                   D1 * ( shellSize - r0 - radius1 ) ) / den1

            a_r_1 = self.D_tot * ( sqrtD_geom * r0 + sqrtD_tot * 
                                   ( shellSize - radius1 ) ) / den1

            assert a_R_1 + a_r_1 * D1_factor + radius1 >= \
                a_R_1 + a_r_1 * D2_factor + radius2

            assert abs( a_R_1 + a_r_1 * D1_factor + radius1 - shellSize ) \
                < 1e-12 * shellSize

            self.a_r = a_r_1
            self.a_R = a_R_1
        else:
            den2 = qrrtD15D2 + D2 * ( sqrtD_geom + sqrtD_tot )

            a_R_2 = sqrtD_geom * ( D1 * ( shellSize - radius2 ) + 
                                   D2 * ( shellSize - r0 - radius2 ) ) / den2

            a_r_2 = self.D_tot * ( sqrtD_geom * r0 + sqrtD_tot * 
                                   ( shellSize - radius2 ) ) / den2

            assert a_R_2 + a_r_2 * D2_factor + radius2 >= \
                a_R_2 + a_r_2 * D1_factor + radius1

            assert abs( a_R_2 + a_r_2 * D2_factor + radius2 - shellSize ) \
                < 1e-12 * shellSize

            self.a_r = a_r_2
            self.a_R = a_R_2

        #log.debug( 'a %g, r %g, R %g r0 %g' % 
        #           ( shellSize, self.a_r, self.a_R, r0 ) )
        #log.debug( 'tr %g, tR %g' % 
        #           ( ( ( self.a_r - r0 ) / math.sqrt(6 * self.D_tot))**2,\
        #                 (self.a_R / math.sqrt( 6*self.D_geom ))**2 ) )
        assert self.a_r > 0
        assert self.a_r > r0, '%g %g' % ( self.a_r, r0 )
        assert self.a_R > 0 or ( self.a_R == 0 and ( D1 == 0 or D2 == 0 ) )

        ### Either com leaves it's PD.
        # draw t_R
        try:
            self.t_R = self.drawTime_single( self.a_R )
        except Exception, e:
            raise Exception, 'sgf.drawTime() failed; %s; %s' %\
                ( str( e ), self.sgf.dump() )

        ### Or r leaves it's PD.
        # draw t_r
        try:
            self.t_r = self.drawTime_pair( r0, self.a_r )
        except Exception, e:
            raise Exception, \
                'pgf.drawTime() failed; %s; r0=%g, %s' % \
                ( str( e ), r0, self.pgf.dump() )


        ### Or a single reacts.
        # draw t_reaction
        t_reaction1 = self.single1.drawReactionTime()
        t_reaction2 = self.single2.drawReactionTime()

        if t_reaction1 < t_reaction2:
            self.t_single_reaction = t_reaction1
            self.reactingsingle = self.single1
        else:
            self.t_single_reaction = t_reaction2
            self.reactingsingle = self.single2

        ### Choose first event.
        ### Here Pair.dt is set!
        self.dt = min( self.t_R, self.t_r, self.t_single_reaction )

        assert self.dt >= 0
        #log.debug( 'dt %g, t_R %g, t_r %g' % 
        #           ( self.dt, self.t_R, self.t_r ) )

        if self.dt == self.t_r:  # type = 0 (REACTION) or 1 (ESCAPE_r)
            try:
                self.eventType = self.drawEventType( r0, self.t_r, self.a_r )
            except Exception, e:
                raise Exception,\
                    'pgf.drawEventType() failed; %s; r0=%g, %s' %\
                    ( str( e ), r0, self.pgf.dump() )

        elif self.dt == self.t_R: # type = ESCAPE_R (2)
            self.eventType = 2
        elif self.dt == self.t_single_reaction:  # type = single reaction (3)
            self.eventType = 3 
        else:
            raise NeverGetHere

        #assert False


    ### Called by determineNextEvent().
    def drawTime_single( self, a ):
        self.sgf.seta( a )
        rnd = numpy.random.uniform()
        return self.sgf.drawTime( rnd )


    ### Called by determineNextEvent().
    def drawTime_pair( self, r0, a ):
        self.pgf.seta( a )
        rnd = numpy.random.uniform()
        #print 'r0 = ', r0, ', rnd = ', rnd[1],\
        #    self.pgf.dump()
        return self.pgf.drawTime( rnd, r0 )


    ### Called by determineNextEvent().
    def drawEventType( self, r0, t, a ):
        rnd = numpy.random.uniform()
        self.pgf.seta( a )
        ### Todo, this C++ function.
        return self.pgf.drawEventType( rnd, r0, t )


    ### Called by firePair() and breakUpPair().
    ### Draw R.
    def drawR_single( self, t, a ):

        self.sgf.seta( a )

        rnd = numpy.random.uniform()
        try:
            r = self.sgf.drawR( rnd, t )
            while r > self.a_R: # redraw; shouldn't happen often
                # Todo. Why could this ever happen? Because we didn't correct
                # for radius of particle (like getMobilitySize for single's)?
                log.info( 'drawR_single: redraw' )
                rnd = numpy.random.uniform()
                r = self.sgf.drawR( rnd, t )
        except Exception, e:
            # BUG. rnd is not subscribtable.
            raise Exception,\
                'gf.drawR_single() failed; %s; rnd= %g, t= %g, %s' %\
                ( str( e ), rnd[2], t, self.sgf.dump() )

        return r


    '''
    Draw r for the pair inter-particle vector.
    '''
    ### Called by firePair() and breakUpPair().
    ### Draw r.
    ### What's the point of passing those a's around, which always is self.a_r, 
    ### if we use self.a_r directly in for example choosePairGreensFunction()
    ### anyway?
    def drawR_pair( self, r0, t, a ):

        gf = self.choosePairGreensFunction( r0, t )

        if hasattr( gf, 'seta' ):  # FIXME: not clean
            gf.seta( a )

        rnd = numpy.random.uniform()
        try:
            r = gf.drawR( rnd, r0, t )
            # redraw; shouldn't happen often
            while r >= self.a_r or r <= self.sigma: 
                # Todo. Why could this ever happen? Because we didn't correct
                # for radius of particle (like getMobilitySize for single's)?
                # Is it because we use approximate greensfunctions?
                log.info( 'drawR_pair: redraw' )
                #self.sim.rejectedMoves += 1  #FIXME:
                rnd = numpy.random.uniform()
                r = gf.drawR( rnd, r0, t )
        except Exception, e:
            raise Exception,\
                'gf.drawR_pair() failed; %s; rnd= %g, r0= %g, t= %g, %s' %\
                ( str( e ), rnd, r0, t, gf.dump() )


        return r


    '''
    Draw theta for the pair inter-particle vector.
    '''
    ### Called by firePair() and breakUpPair().
    def drawTheta_pair( self, rnd, r, r0, t, a ):

        #print 'r ', r, 'r0 ', r0, 't ', t, 'a ', a
        gf = self.choosePairGreensFunction( r0, t )

        #print gf
        if hasattr( gf, 'seta' ):  # FIXME: not clean
            gf.seta( a )

        try:
            theta = gf.drawTheta( rnd, r, r0, t )
        except Exception, e:
            raise Exception,\
                'gf.drawTheta() failed; %s; rnd= %g, r= %g, r0= %g, t=%g, %s' %\
                ( str( e ), rnd, r, r0, t, gf.dump() )

        return theta


    ### Called by firePair() and breakUpPair().
    def checkNewpos( self, pos1, pos2, com ):

        species1 = self.single1.particle.species
        species2 = self.single2.particle.species

        oldCoM = com
        
        # debug: check if the new positions are valid:
        newDistance = distance_Simple( pos1, pos2 )
        particleRadius12 = species1.radius + species2.radius

        # check 1: particles don't overlap.
        if newDistance <= particleRadius12:
            log.info( 'rejected move: radii %g, particle distance %g',
                          ( species1.radius + species2.radius, newDistance ) )
            log.debug( 'DEBUG: dt %g, pos1 %s, pos2 %s' %
                           ( self.dt, str( pos1 ), str( pos2 ) ) )
            raise RuntimeError, 'New particles overlap'

        # check 2: particles within mobility radius.
        d1 = self.distance( oldCoM, pos1 ) + species1.radius
        d2 = self.distance( oldCoM, pos2 ) + species2.radius
        if d1 > self.radius or d2 > self.radius:
            raise RuntimeError, \
                'New particle(s) out of protective sphere. %s' % \
                'radius = %g, d1 = %g, d2 = %g ' % ( self.radius, d1, d2 )
                
        

        return True


    def check( self ):
        pass

    def __str__( self ):
        buf = 'Pair( ' + str(self.single1.particle) +\
              ', ' + str(self.single2.particle) + ' )'

        return buf
