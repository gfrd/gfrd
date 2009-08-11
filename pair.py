from utils import *
from _gfrd import *
from shape import *
from domain import *
from log import *


'''
There are 3 types of pairs:
    * SphericalPair
    * PlanarSurfacePair
    * CylindricalSurfacePair

Todo. Add more comments.
'''
class Pair( object ):
    '''
    Todo. How does this work exactly?
    CUTOFF_FACTOR is a threshold to choose between the real and approximate
    Green's functions.
    H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7,
    5.6: ~1e-8, 6.0: ~1e-9
    '''
    CUTOFF_FACTOR = 5.6

    '''
    The distance function is passed to be able to take the boundary 
    conditions into account.
    '''
    def __init__( self, single1, single2, shellSize, rt, distFunc, worldSize ):
        self.multiplicity = 2

        '''
        Individual singles are preserved, but removed from shellMatrix and 
        eventScheduler.
        Order single1 and single2 so that D1 < D2.
        '''
        if single1.particle.species.D <= single2.particle.species.D:
            self.single1, self.single2 = single1, single2 
        else:
            self.single1, self.single2 = single2, single1 
        self.singles = [ self.single1, self.single2 ]

        assert single1.surface == single2.surface
        self.surface = single1.surface

        self.rt = rt

        self.distance = distFunc
        self.worldSize = worldSize
        
        particle1 = self.single1.particle
        particle2 = self.single2.particle

        self.biggestParticleRadius = max( particle1.species.radius, particle2.species.radius )

        self.D1, self.D2 = particle1.species.D, particle2.species.D
        # Todo. Is this also correct for 1D and 2D?
        self.D_tot = self.D1 + self.D2
        self.D_geom = math.sqrt( self.D1 * self.D2 )  # geometric mean

        self.eventID = None
        self.lastTime = 0.0
        self.dt = 0.0
        self.eventType = None

        self.sigma = particle1.species.radius + particle2.species.radius
        self.shellSize = shellSize
        # getCoMandIV calls cyclicTranspose on particle2.pos, so both CoM as 
        # well as IV are correct (periodic boundary conditions taken into 
        # account). 
        self.CoM, self.IV = self.getCoMandIV()
        self.pairDistance = distFunc( particle1.pos, particle2.pos )
        # Todo: why discrepancy of about 1e-15 here?
        assert abs(self.pairDistance - length( self.IV )) < 1e-10, \
                    '%.15g != %.15g' % (self.pairDistance, length( self.IV ))


    def __del__( self ):
        #log.debug( 'del %s' % str( self ) )
        pass


    def initialize( self, t ):
        self.lastTime = t
        # Don't reset radius. Pairs don't move.
        self.dt = 0
        self.eventType = None


    def setPos( self, pos ):
        raise Exception
        self.shellList[0].origin = pos
    def getPos( self ):
        raise Exception
        return self.shellList[0].origin
    pos = property( getPos, setPos )


    def getD( self ):
        return self.D_tot #FIXME: is this correct?


    '''
    This method returns the radius from its CoM that this Pair must reserve
    to remain mobile.
    '''
    def getMinRadius( self ):
        minRadius = max( pairDistance * self.D1 /
                         self.D_tot + self.single1.getMinRadius(),
                         pairDistance * self.D2 /
                         self.D_tot + self.single2.getMinRadius() )
        return minRadius


    '''
    Calculate and return the "Center of Mass" (== CoM) of this pair.
    '''
    def getCoMandIV( self ):
        pos1 = self.single1.particle.pos
        pos2 = self.single2.particle.pos

        # Todo. Why the fixme?
        pos2t = cyclicTranspose( pos2, pos1, self.worldSize ) #FIXME
        CoM = ( pos1 * self.D2 + pos2t * self.D1 ) / self.D_tot
        IV = pos2t - pos1
        return CoM % self.worldSize, IV


    '''
    Determine a_r and a_R.

    Todo. Check all this and make sure the radius of the largest particle is 
    subtracted from a_r and a_R to get values of the same type as 
    getMobilityRadius() for singles. 
    
          Make dimension specific.
    '''
    def determineRadii( self ):
        single1 = self.single1
        single2 = self.single2
        radius1 = single1.particle.species.radius
        radius2 = single2.particle.species.radius

        D1 = self.D1
        D2 = self.D2

        D1_factor = D1 / self.D_tot
        D2_factor = D2 / self.D_tot

        shellSize = self.shellSize / SAFETY  # FIXME:

        sqrtD_tot = math.sqrt( self.D_tot )
        sqrtD_geom = math.sqrt( self.D_geom )

        pairDistance = self.pairDistance

        assert pairDistance >= self.sigma, \
            '%s;  pairDistance %g < sigma %g' % ( self, pairDistance, self.sigma )

        # equalize expected mean t_r and t_R.

        qrrtD1D25 = ( D1    * D2**5 ) ** 0.25
        qrrtD15D2 = ( D1**5 * D2 ) ** 0.25

        if qrrtD15D2 * pairDistance + ( qrrtD15D2 + qrrtD1D25 ) * radius1 \
                + D1 * ( sqrtD_tot * ( shellSize - radius2 ) 
                         - sqrtD_geom * radius2 )\
                - D2 * ( sqrtD_geom * pairDistance + sqrtD_tot * 
                         ( shellSize - radius1 ) )\
                         - qrrtD1D25 * radius2 >= 0:

            den1 = qrrtD1D25 + D1 * ( sqrtD_geom + sqrtD_tot )

            a_R_1 = sqrtD_geom * ( D2 * ( shellSize - radius1) + 
                                   D1 * ( shellSize - pairDistance - radius1 ) ) / den1

            a_r_1 = self.D_tot * ( sqrtD_geom * pairDistance + sqrtD_tot * 
                                   ( shellSize - radius1 ) ) / den1

            assert a_R_1 + a_r_1 * D1_factor + radius1 >= \
                a_R_1 + a_r_1 * D2_factor + radius2

            assert abs( a_R_1 + a_r_1 * D1_factor + radius1 - shellSize ) \
                < 1e-12 * shellSize

            a_r = a_r_1
            a_R = a_R_1
        else:
            den2 = qrrtD15D2 + D2 * ( sqrtD_geom + sqrtD_tot )

            a_R_2 = sqrtD_geom * ( D1 * ( shellSize - radius2 ) + 
                                   D2 * ( shellSize - pairDistance - radius2 ) ) / den2

            a_r_2 = self.D_tot * ( sqrtD_geom * pairDistance + sqrtD_tot * 
                                   ( shellSize - radius2 ) ) / den2

            assert a_R_2 + a_r_2 * D2_factor + radius2 >= \
                a_R_2 + a_r_2 * D1_factor + radius1

            assert abs( a_R_2 + a_r_2 * D2_factor + radius2 - shellSize ) \
                < 1e-12 * shellSize

            a_r = a_r_2
            a_R = a_R_2

        #log.debug( 'a %g, r %g, R %g pairDistance %g' % 
        #           ( shellSize, a_r, a_R, pairDistance ) )
        #log.debug( 'tr %g, tR %g' % 
        #           ( ( ( a_r - pairDistance ) / math.sqrt(6 * self.D_tot))**2,\
        #                 (a_R / math.sqrt( 6*self.D_geom ))**2 ) )
        assert a_r > 0
        assert a_r > pairDistance, '%g %g' % ( a_r, pairDistance )
        assert a_R > 0 or ( a_R == 0 and ( D1 == 0 or D2 == 0 ) )
        return a_R, a_r

     
    # Returns a (dt, activeDomain) tuple.
    def drawEscapeOrReactionTime( self ):
        '''
        Note: we are not deciding yet if this is an escape or a pair 
        reaction, since we aren't calling domain.drawEventType() yet. We 
        postpone it to the very last minute (when this event is executed 
        in firePair), and memorize the activeDomain (CoM or IV) like this.

        Also note that in case the dt for single reaction is actually smaller 
        below, and this single will get a single reaction event, 
        self.activeDomain won't be used at all, as long as you make sure 
        reaction events are taken care of before escape events in firePair.

        Last note: 
        * EventType.REACTION means *single* reaction,
        * EventType.ESCAPE can still mean any of CoM or IV escape or *pair* 
          reaction.
        '''
        return min( (d.drawTime(), d) for d in self.domains )


    # Returns a (dt, single) tuple.
    def drawSingleReactionTime( self ):
        return min( ((single.drawReactionTime(), single) for single in self.singles ))


    def checkNewpos( self, pos1, pos2 ):
        species1 = self.single1.particle.species
        species2 = self.single2.particle.species

        oldCoM = self.CoM
        
        # debug: check if the new positions are valid:
        newDistance = distance_Simple( pos1, pos2 )
        particleRadius12 = species1.radius + species2.radius

        # check 1: particles don't overlap.
        if newDistance <= particleRadius12:
            raise Stop( 'New particles overlap' )

        # check 2: particles within mobility radius.
        d1 = self.distance( oldCoM, pos1 ) + species1.radius
        d2 = self.distance( oldCoM, pos2 ) + species2.radius
        if d1 > self.shellSize or d2 > self.shellSize:
            raise Stop( 'New particle(s) out of protective sphere. %s' % \
                'radius = %g, d1 = %g, d2 = %g ' % ( self.shellSize, d1, d2 ) )
        return True


    def check( self ):
        pass


    def __str__( self ):
        buf = '( ' + str(self.single1.particle) +\
              ', ' + str(self.single2.particle) + ' )'
        return buf


##############################################################################
'''
2 Particles inside a (spherical) shell not on any surface.
'''
class SphericalPair( Pair ):
    def __init__( self, single1, single2, shellSize, rt, distFunc, worldSize ):
        Pair.__init__( self, single1, single2, shellSize, rt, distFunc, worldSize )

        # Set shellSize directly, without getMinRadius() step. Don't set 
        # radius again from initialize(). Pairs don't move.
        self.shellList = [ Sphere( self.CoM, shellSize ) ]

        a_R, self.a_r = self.determineRadii()

        # Green's function for centre of mass inside absorbing sphere.
        sgf = FirstPassageGreensFunction( self.D_geom )
        comDomain = RadialDomain( a_R, sgf )

        # Green's function for interparticle vector inside absorbing sphere.  
        # This exact solution is used for drawing times.
        self.pgf = FirstPassagePairGreensFunction( self.D_tot, self.rt.k, 
                                                        self.sigma )
        ivDomain = CompositeDomain( self.sigma, self.pairDistance, self.a_r, self.pgf )

        self.domains = [ comDomain, ivDomain ]


    '''
    Calculate new positions of the pair particles using a new center-of-mass, 
    a new inter-particle vector, and an old inter-particle vector.
    '''
    def drawNewPositions( self, dt ):
        newCoM, newIV = self.drawNewCoM( dt ), self.drawNewIV( dt )
        #FIXME: need better handling of angles near zero and pi.

        '''
        Rotate the new interparticle vector along the rotation axis that is 
        perpendicular to both the z-axis and the original interparticle vector 
        for the angle between these.
        
        the rotation axis is a normalized cross product of the z-axis and the 
        original vector.
        rotationAxis = crossproduct( [ 0,0,1 ], IV )
        '''

        oldIV = self.IV
        angle = vectorAngleAgainstZAxis( oldIV )
        if angle % numpy.pi != 0.0:
            rotationAxis = crossproductAgainstZAxis( oldIV )
            rotationAxis = normalize( rotationAxis )
            rotated = rotateVector( newIV,
                                    rotationAxis,
                                    angle )
        elif angle == 0.0:
            rotated = newIV
        else:
            rotated = numpy.array( [ newIV[0], newIV[1], - newIV[2] ] )
            #rotated = newIV * numpy.array( [ 1, 1, -1 ] )

        newpos1 = newCoM - rotated * ( self.D1 / self.D_tot )
        newpos2 = newCoM + rotated * ( self.D2 / self.D_tot )
        return newpos1, newpos2


    def drawNewCoM( self, dt ):
        r_R = self.domains[0].drawPosition( dt )
        return self.CoM + randomVector( r_R )


    def drawNewIV( self, dt ):  
        gf = self.choosePairGreensFunction( dt )
        r, theta = self.domains[1].drawPosition( gf, dt )
        newInterParticleS = numpy.array([r, theta, numpy.random.random()*2*Pi])
        return sphericalToCartesian( newInterParticleS )


    '''
    Decide if we need to use an approximate Green's function for drawing 
    positions. Really needed for small dt for example.
    ro is interparticle distance.
    '''
    def choosePairGreensFunction( self, dt ):
        distanceFromSigma = self.pairDistance - self.sigma
        distanceFromShell = self.a_r - self.pairDistance

        thresholdDistance = Pair.CUTOFF_FACTOR * \
            math.sqrt( 6.0 * self.D_tot * dt )

        if distanceFromSigma < thresholdDistance:
            if distanceFromShell < thresholdDistance:
                # near both a and sigma;
                # use FirstPassagePairGreensFunction
                #log.debug( 'GF: normal' )
                pgf = self.pgf
                pgf.seta( self.a_r ) # Don't forget.
            else:
                # near sigma; use BasicPairGreensFunction (has no a).
                #log.debug( 'GF: only sigma' )
                pgf = BasicPairGreensFunction( self.D_tot, self.rt.k, 
                                                self.sigma )
        else:
            if distanceFromShell < thresholdDistance:
                # near a;
                #log.debug( 'GF: only a' )
                pgf = FirstPassageNoCollisionPairGreensFunction( self.D_tot )
                pgf.seta( self.a_r ) # Don't forget.
                
            else:
                # distant from both a and sigma;
                # use FreePairGreensFunction (has no a).
                #log.debug( 'GF: free' )
                pgf = FreePairGreensFunction( self.D_tot )

        return pgf


    def __str__( self ):
        return 'SphericalPair' + Pair.__str__( self )


'''
2 Particles inside a (cylindrical) shell on a PlanarSurface. (Hockey pucks).
'''
class PlanarSurfacePair( Pair ):
    def __init__( self, single1, single2, shellSize, rt, distFunc, worldSize ):
        Pair.__init__( self, single1, single2, shellSize, rt, distFunc, worldSize )

        '''
        Hockey pucks stick a bit out of the surface, so that if the 
        biggest particle would undergo an unbinding reaction, it can still be 
        placed within the hockey puck (and thus it won't interfere with other 
        shells.

        Set radius=shellSize directly, without getMinRadius() step. Don't set 
        radius again from initialize(). Pairs don't move.
        '''
        self.shellList = [ Cylinder( self.CoM, shellSize, self.surface.unitZ, self.surface.Lz * SAFETY + 2*self.biggestParticleRadius ) ]

        a_R, self.a_r = self.determineRadii()

        # Green's function for centre of mass inside absorbing sphere.
        sgf = FirstPassageGreensFunction( self.D_geom )
        comDomain = RadialDomain( a_R, sgf )

        # Green's function for interparticle vector inside absorbing sphere.  
        # This exact solution is used for drawing times.
        self.pgf = FirstPassagePairGreensFunction2D( self.D_tot, self.rt.k, 
                                                        self.sigma )
        ivDomain = CompositeDomain( self.sigma, self.pairDistance, self.a_r, self.pgf )

        self.domains = [ comDomain, ivDomain ]


    '''
    Calculate new positions of the pair particles using a new center-of-mass, 
    a new inter-particle vector, and an old inter-particle vector.
    '''
    def drawNewPositions( self, dt ):
        newCoM, (r, theta) = self.drawNewCoM( dt ), self.drawNewIV( dt )
        unitX = self.surface.unitX
        unitY = self.surface.unitY
        angle = vectorAngle( unitX, self.IV )
        # Todo. Test if nothing changes when theta=0.
        newAngle = angle + theta
        
        rotated = r * math.cos(newAngle) * unitX + r * math.sin(newAngle) * unitY

        newpos1 = newCoM - rotated * ( self.D1 / self.D_tot )
        newpos2 = newCoM + rotated * ( self.D2 / self.D_tot )
        return newpos1, newpos2


    def drawNewCoM( self, dt ):
        r_R = self.domains[0].drawPosition( dt )
        x, y = randomVector2D( r_R )
        return self.CoM + x * self.surface.unitX + y * self.surface.unitY


    def drawNewIV( self, dt ):
        r, theta = self.domains[1].drawPosition( self.pgf, dt )
        assert r > self.sigma and r <= self.a_r
        return r, theta


    def __str__( self ):
        return 'PlanarSurfacePair' + Pair.__str__( self )


'''
2 Particles inside a (cylindrical) shell on a CylindricalSurface. (Rods).
'''
class CylindricalSurfacePair( Pair ):
    def __init__( self, single1, single2, shellSize, rt, distFunc, worldSize ):
        Pair.__init__( self, single1, single2, shellSize, rt, distFunc, worldSize )

        # Set shellSize directly, without getMinRadius() step. Don't set 
        # radius again from initialize(). Pairs don't move.
        self.shellList = [ Cylinder( self.CoM, self.surface.radius * SAFETY + 2*self.biggestParticleRadius, self.surface.unitZ, shellSize ) ]

        a_R, self.a_r = self.determineRadii()

        # Green's function for centre of mass inside absorbing sphere.
        sgf = FirstPassageGreensFunction1D( self.D_geom )
        # a_R can be used as is for cartesian domain.
        comDomain = CartesianDomain( 0, a_R, sgf )

        # Green's function for interparticle vector.
        self.pgf = FirstPassageGreensFunction1DRad( self.D_tot, self.rt.k )
        # Calculate a and r0 for a cartesian domain. A bit tricky. Needed 
        # because normally iv is used with an r-theta domain, and then r goes 
        # from sigma (radiating boundary) to a_r (absorbing boundary). Now, 
        # cartesian domain goes from -a_cartesian (radiating boundary) to 
        # a_cartesian (absorbing boundary).
        a_cartesian = (self.a_r - self.sigma ) / 2
        r0_cartesian = self.pairDistance - self.sigma - a_cartesian
        ivDomain = CartesianDomain( r0_cartesian, a_cartesian, self.pgf )

        self.domains = [ comDomain, ivDomain ]


    def drawNewPositions( self, dt ):
        newCoM, newIV = self.drawNewCoM( dt ), self.drawNewIV( dt )
        newpos1 = newCoM - newIV * ( self.D1 / self.D_tot )
        newpos2 = newCoM + newIV * ( self.D2 / self.D_tot )
        return newpos1, newpos2


    def drawNewCoM( self, dt ):
        r_R = self.domains[0].drawPosition( dt )
        return self.CoM + r_R * self.surface.unitZ


    def drawNewIV( self, dt ):
        r_cartesian = self.domains[1].drawPosition( dt )
        # Convert back from cartesian domain (-a_cartesian to a_cartesian) to 
        # something that can be used for the length of the iv vector (sigma to 
        # a_r.
        a_cartesian = self.domains[1].a
        r = r_cartesian + self.sigma + a_cartesian
        assert r > self.sigma and r <= self.a_r
        # Note: using self.surface.unitZ here might accidently interchange the 
        # particles.
        return r * normalize( self.IV )


    def __str__( self ):
        return 'CylindricalSurfacePair' + Pair.__str__( self )


