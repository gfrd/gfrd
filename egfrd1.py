#!/usr/env python

import weakref
import math
import numpy

from gfrdbase import *

from single import *
from pair import *
from multi import *
from shape import Sphere, Cylinder

class Delegate( object ):

    def __init__( self, obj, method ):
        self.obj = weakref.proxy( obj )
        self.method = method

    def __call__( self, *arg ):
        return self.method( self.obj, *arg )

# Notes:
# EGFRDSimulator.distance takes into account the periodic boundary conditions.  
# Very handy!

class EGFRDSimulator1( ParticleSimulatorBase ):
    
    def __init__( self ):

        self.sphereMatrix = SphereMatrix()
        self.cylinderMatrix = CylinderMatrix()
        #self.sm2 = pObjectMatrix()

        ParticleSimulatorBase.__init__( self )

        self.MULTI_SHELL_FACTOR = 0.05
        self.SINGLE_SHELL_FACTOR = 0.1

        self.isDirty = True
        self.scheduler = EventScheduler()

        self.smallT = 1e-8  # FIXME: is this ok?

        self.userMaxShellSize = INF

        self.reset()


    def setWorldSize( self, size ):
        if isinstance( size, list ) or isinstance( size, tuple ):
            size = numpy.array( size )

        ParticleSimulatorBase.setWorldSize( self, size )
        self.sphereMatrix.setWorldSize( size )
        self.cylinderMatrix.setWorldSize( size )


    def setMatrixSize( self, size ):
        ParticleSimulatorBase.setMatrixSize( self, size )
        self.sphereMatrix.setMatrixSize( size )
        self.cylinderMatrix.setMatrixSize( size )


    def getMatrixCellSize( self ):
        assert self.sphereMatrix.cellSize == self.cylinderMatrix.cellSize
        return self.sphereMatrix.cellSize


    def getNextTime( self ):

        if self.scheduler.getSize() == 0:
            return self.t

        return self.scheduler.getTopTime()

    def setUserMaxShellSize( self, size ):

        self.userMaxShellSize = size

    def getUserMaxShellSize( self ):

        return self.userMaxShellSize

    def getMaxShellSize( self ):

        return min( self.getMatrixCellSize() * .5 / SAFETY,
                    self.userMaxShellSize )

    def reset( self ):
        self.t = 0.0
        self.dt = 0.0
        self.stepCounter = 0
        self.zeroSteps = 0
        self.rejectedMoves = 0
        self.reactionEvents = 0
        self.lastEvent = None
        self.lastReaction = None

        self.isDirty = True
        #self.initialize()
        

    def initialize( self ):
        ParticleSimulatorBase.initialize( self )

        self.setAllRepulsive()

        self.scheduler.clear()
        self.sphereMatrix.clear()
        self.cylinderMatrix.clear()

        # Get particle data from particleMatrix, because we need to know the 
        # surface they are on.
        particles, _ = self.particleMatrix.getNeighbors( numpy.array([0,0,0]) )
        for p in particles:
            # Get reference to real Particle. I think this is needed.
            particle = Particle( p.species, serial=p.serial, surface=p.surface )
            single = self.createSingle( particle )
            self.addToShellMatrix( single )
            self.addSingleEvent( single )

        self.isDirty = False


    def stop( self, t ):

        log.info( 'stop at %g' % t )

        if self.t == t:
            return

        if t >= self.scheduler.getTopEvent().getTime():
            raise RuntimeError, 'Stop time >= next event time.'

        if t < self.t:
            raise RuntimeError, 'Stop time < current time.'

        self.t = t
        
        scheduler = self.scheduler
        
        nonSingleList = []

        # first burst all Singles.
        for i in range( scheduler.getSize() ):
            obj = scheduler.getEventByIndex(i).getArg()
            if isinstance( obj, Pair ) or isinstance( obj, Multi ):
                nonSingleList.append( obj )
            elif isinstance( obj, Single ):
                log.debug( 'burst %s, lastTime= %g' % 
                           ( str( obj ), obj.lastTime ) )
                self.burstSingle( obj )
            else:
                assert False, 'do not reach here'


        # then burst all Pairs and Multis.
        log.debug( 'burst %s' % nonSingleList )
        self.burstObjs( nonSingleList )

        self.dt = 0.0


    def step( self ):

        self.lastReaction = None

        if self.isDirty:
            self.initialize()
            
        #if self.stepCounter % 100 == 0:
        #    self.check()
        
        self.stepCounter += 1

        event = self.scheduler.getTopEvent()
        self.t, self.lastEvent = event.getTime(), event.getArg()

        log.info( '\n%d: t=%g dt=%g\nevent=%s reactions=%d rejectedmoves=%d' 
                      % ( self.stepCounter, self.t, self.dt, self.lastEvent, 
                          self.reactionEvents, self.rejectedMoves ) )
        
        self.scheduler.step()

        nextTime = self.scheduler.getTopTime()
        self.dt = nextTime - self.t


        # assert if not too many successive dt=0 steps occur.
        if self.dt == 0:
            self.zeroSteps += 1
            if self.zeroSteps >= max( self.scheduler.getSize() * 3, 10 ):
                raise RuntimeError, 'too many dt=zero steps.  simulator halted?'
        else:
            self.zeroSteps = 0


        assert self.scheduler.getSize() != 0


    def createMulti( self ):

        multi = Multi( self )

        #multi.initialize( self.t )

        return multi


    def moveSingle( self, single, pos ):
        single.pos = pos
        self.updateOnParticleMatrix( single.particle, pos )


    def addToShellMatrix( self, obj ):
        for i, shell in enumerate( obj.shellList ):
            key = (obj, i)
            if isinstance( shell, Sphere ):
                self.sphereMatrix.add( key, shell.origin, shell.size )
            elif isinstance( shell, Cylinder ):
                self.cylinderMatrix.add( key, shell )
            else: raise KeyError, 'Objecttype does not exit'


    def removeFromShellMatrix( self, obj ):
        for i, shell in enumerate( obj.shellList ):
            key = (obj, i)
            # Check for type of shell -> remove key (is that the way to do 
            # it?)
            if isinstance( shell, Sphere ):
                self.sphereMatrix.remove( key )
            elif isinstance( shell, Cylinder ):
                self.cylinderMatrix.remove( key )
            else: raise KeyError, 'Objecttype does not exit'


    def updateShellMatrix( self, obj ):
        for i, shell in enumerate( obj.shellList ):
            key = (obj, i)
            if isinstance( shell, Sphere ):
                self.sphereMatrix.update( key, shell.origin, shell.size )
            elif isinstance( shell, Cylinder ):
                self.cylinderMatrix.update( key, shell )
            else: raise KeyError, 'Objecttype does not exit'


    def addEvent( self, t, func, arg ):

        return self.scheduler.addEvent( t, func, arg )


    def removeEvent( self, event ):

        self.scheduler.removeEvent( event.eventID )


    def updateEvent( self, t, event ):
        self.scheduler.updateEventTime( event.eventID, t )


    def burstObj( self, obj ):
        
        log.info( 'bursting %s' % str( obj ) )

        if isinstance( obj, Single ):
            self.burstSingle( obj )
            return [obj,]
        elif isinstance( obj, Pair ):  # Pair
            single1, single2 = self.burstPair( obj )
            self.removeEvent( obj )
            self.addSingleEvent( single1 )
            self.addSingleEvent( single2 )
            return [ single1, single2 ]
        else:  # Multi
            bursted = self.burstMulti( obj )
            self.removeEvent( obj )
            return bursted


    def burstObjs( self, objs ):

        bursted = []
        for obj in objs:
            b = self.burstObj( obj )
            bursted.extend( b )

        return bursted

    def clearVolume( self, pos, radius, ignore=[] ):

        neighbors = self.getNeighborsWithinRadiusNoSort( pos, radius, ignore )

        return self.burstObjs( neighbors )


    def burstNonMultis( self, neighbors ):

        bursted = []

        for obj in neighbors:
            if not isinstance( obj, Multi ):
                b = self.burstObj( obj )
                bursted.extend( b )
            else:
                bursted.append( obj )

        return bursted



    def restoreSingleShells( self, singles ):

        for single in singles:
            assert single.isReset()
            c, d = self.getClosestObj( single.pos, ignore = [single,] )

            self.updateSingle( single, c, d )
            self.updateEvent( self.t + single.dt, single )
            log.debug( 'restore shell %s %g dt %g closest %s %g' %
                       ( single, single.size, single.dt, c, d ) )


    def calculateSingleShellSize( self, single, closest, 
                                  distance, shellDistance ):

        minSize1 = single.getMinSize()
        D1 = single.getD()

        if D1 == 0:
            return minSize1

        D2 = closest.getD()
        minSize2 = closest.getMinSize()
        minSize12 = minSize1 + minSize2
        sqrtD1 = math.sqrt( D1 )
            
        shellSize = min( sqrtD1 / ( sqrtD1 + math.sqrt( D2 ) )
                         * ( distance - minSize12 ) + minSize1,
                         shellDistance )
        shellSize /= SAFETY
        shellSize = max( shellSize, minSize1 ) # not smaller than the radius

        return shellSize


    def updateSingle( self, single, closest, distanceToShell ): 

        if isinstance( closest, Single ):
            distanceToClosest = self.distance( single.pos, closest.pos )
            shellSize = self.calculateSingleShellSize( single, closest, 
                                                       distanceToClosest,
                                                       distanceToShell )
        else:  # Pair or Multi
            shellSize = distanceToShell / SAFETY
            shellSize = max( shellSize, single.getMinSize() )

        shellSize = min( shellSize, self.getMaxShellSize() )

        single.setSize( shellSize )
        single.determineNextEvent( self.t )
        self.updateShellMatrix( single )



    def firePair( self, pair ):

        assert self.checkObj( pair )

        log.info( 'fire: %s eventType %s' % ( pair, pair.eventType ) )

        particle1 = pair.single1.particle
        particle2 = pair.single2.particle
        
        oldInterParticle = particle2.pos - particle1.pos
        oldCoM = pair.getCoM()
        self.applyBoundary( oldCoM )

        # Three cases:
        #  0. Reaction
        #  1. Escaping through a_r.
        #  2. Escaping through a_R.
        #  3. Single reaction 

        # First handle single reaction case.
        if pair.eventType == 3:

            reactingsingle = pair.reactingsingle

            log.info( 'pair: single reaction %s' % str( reactingsingle ) )

            if reactingsingle == pair.single1:
                theothersingle = pair.single2
            else:
                theothersingle = pair.single1

            self.burstPair( pair )

            self.addSingleEvent( theothersingle )

            try:
                self.removeFromShellMatrix( reactingsingle )
                self.fireSingleReaction( reactingsingle )
            except NoSpace:
                self.addToShellMatrix( reactingsingle )
                self.rejectedMoves += 1
                reactingsingle.dt = 0
                self.addSingleEvent( reactingsingle )

            pair.dt = -INF
            return pair.dt
        


        #
        # 0. Reaction
        #
        if pair.eventType == EventType.REACTION:

            log.info( 'reaction' )

            if len( pair.rt.products ) == 1:
                
                species3 = pair.rt.products[0]

                rnd = numpy.random.uniform( size=2 )

                # calculate new R
            
                r_R = pair.drawR_single( pair.dt, pair.a_R )
            
                displacement_R_S = [ r_R, rnd[0] * Pi, rnd[1] * 2 * Pi ]
                displacement_R = sphericalToCartesian( displacement_R_S )
                newCoM = oldCoM + displacement_R
                
                assert self.distance( oldCoM, newCoM ) + species3.radius <\
                    pair.radius

                #FIXME: SURFACE
                currentSurface = particle1.surface
                self.applyBoundary( newCoM )

                self.removeParticle( particle1 )
                self.removeParticle( particle2 )

                particle = self.createParticle( species3, newCoM, currentSurface  )
                newsingle = self.createSingle( particle )
                self.addToShellMatrix( newsingle )
                self.addSingleEvent( newsingle )

                self.reactionEvents += 1

                self.lastReaction = Reaction( pair.rt, [particle1, particle2],
                                              [particle] )

                log.info( 'product; %s' % str( newsingle ) )

            else:
                raise NotImplementedError,\
                      'num products >= 2 not supported.'

            self.removeFromShellMatrix( pair )

            pair.dt = -INF
            return pair.dt


        #
        # Escape 
        #

        r0 = self.distance( particle1.pos, particle2.pos )

        # 1 Escaping through a_r.
        if pair.eventType == EventType.ESCAPE:

            log.debug( 'r0 = %g, dt = %g, %s' %
                           ( r0, pair.dt, pair.pgf.dump() ) )
            
            rnd = numpy.random.uniform( size=4 )

            # calculate new R
            
            r_R = pair.drawR_single( pair.dt, pair.a_R )
                
            displacement_R_S = [ r_R, rnd[0] * Pi, rnd[1] * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = oldCoM + displacement_R

            # calculate new r
            theta_r = pair.drawTheta_pair( rnd[2], pair.a_r, r0, pair.dt, 
                                           pair.a_r )
            phi_r = rnd[3] * 2 * Pi
            newInterParticleS = numpy.array( [ pair.a_r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )
                
            newpos1, newpos2 = pair.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )
            self.applyBoundary( newpos1 )
            self.applyBoundary( newpos2 )


        # 2 escaping through a_R.
        elif pair.eventType == 2:

            rnd = numpy.random.uniform( size = 4 )

            # calculate new r
            log.debug( 'r0 = %g, dt = %g, %s' %
                           ( r0, pair.dt, pair.pgf.dump() ) )
            r = pair.drawR_pair( r0, pair.dt, pair.a_r )
            log.debug( 'new r = %g' % r )
            #assert r >= pair.sigma
            
            theta_r = pair.drawTheta_pair( rnd[0], r, r0, pair.dt, pair.a_r )
            phi_r = rnd[1] * 2*Pi
            newInterParticleS = numpy.array( [ r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )
                
            # calculate new R
            displacement_R_S = [ pair.a_R, rnd[2] * Pi, rnd[3] * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            
            newCoM = oldCoM + displacement_R
                
            newpos1, newpos2 = pair.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )
            self.applyBoundary( newpos1 )
            self.applyBoundary( newpos2 )

        else:
            raise SystemError, 'Bug: invalid eventType.'

        # this has to be done before the following clearVolume()

        self.removeFromShellMatrix( pair )

        assert pair.checkNewpos( newpos1, newpos2, oldCoM )
        assert self.checkOverlap( newpos1, particle1.species.radius,
                                  ignore = [ particle1, particle2 ] )
        assert self.checkOverlap( newpos2, particle2.species.radius,
                                  ignore = [ particle1, particle2 ] )

        single1, single2 = pair.single1, pair.single2

        self.moveSingle( single1, newpos1 )
        self.moveSingle( single2, newpos2 )

        single1.initialize( self.t )
        single2.initialize( self.t )
            
        self.addSingleEvent( single1 )
        self.addSingleEvent( single2 )

        self.addToShellMatrix( single1 )
        self.addToShellMatrix( single2 )

        assert self.checkObj( single1 )
        assert self.checkObj( single2 )

        pair.dt = -INF
        return pair.dt


    def fireMulti( self, multi ):
        
        sim = multi.sim

        sim.step()
        #sim.sync()

        if sim.lastReaction:
            log.info( 'bd reaction' )

            self.breakUpMulti( multi )
            self.reactionEvents += 1
            self.lastReaction = sim.lastReaction
            return -INF

        if sim.escaped:
            log.info( 'multi particle escaped.' )

            self.breakUpMulti( multi )
            return -INF

        #log.info( 'multi stepped %d steps, duration %g, dt = %g' %
        #          ( additionalSteps + 1, sim.t - startT + sim.dt, dt ) )

        return multi.dt


    def breakUpMulti( self, multi ):

        self.removeFromShellMatrix( multi )

        singles = []
        for particle in multi.sim.particleList:
            single = self.createSingle( particle )
            self.addToShellMatrix( single )
            self.addSingleEvent( single )
            singles.append( single )

        return singles


    def burstMulti( self, multi ):
        
        #multi.sim.sync()
        singles = self.breakUpMulti( multi )

        return singles


    def breakUpPair( self, pair ):

        assert self.t >= pair.lastTime
        assert self.t <= pair.lastTime + pair.dt

        dt = self.t - pair.lastTime 

        if dt > 0.0:

            single1 = pair.single1
            single2 = pair.single2
            particle1 = single1.particle
            particle2 = single2.particle

            oldInterParticle = single2.pos - single1.pos
            oldCoM = pair.getCoM()
            r0 = pair.distance( single1.pos, single2.pos )
            
            rnd = numpy.random.uniform( size = 4 )

            # calculate new CoM
            r_R = pair.drawR_single( dt, pair.a_R )
            
            displacement_R_S = [ r_R, rnd[0] * Pi, rnd[1] * 2 * Pi ]
            displacement_R = sphericalToCartesian( displacement_R_S )
            newCoM = oldCoM + displacement_R
            
            # calculate new interparticle
            r_r = pair.drawR_pair( r0, dt, pair.a_r )
            theta_r = pair.drawTheta_pair( rnd[2], r_r, r0, dt, pair.a_r )
            phi_r = rnd[3] * 2 * Pi
            newInterParticleS = numpy.array( [ r_r, theta_r, phi_r ] )
            newInterParticle = sphericalToCartesian( newInterParticleS )

            newpos1, newpos2 = pair.newPositions( newCoM, newInterParticle,
                                                  oldInterParticle )
            self.applyBoundary( newpos1 )
            self.applyBoundary( newpos2 )
            assert self.checkOverlap( newpos1, particle1.species.radius,
                                      ignore = [ particle1, particle2 ] )
                                      
            assert self.checkOverlap( newpos2, particle2.species.radius,
                                      ignore = [ particle1, particle2 ] )
                                      
            assert pair.checkNewpos( newpos1, newpos2, oldCoM )
            self.moveSingle( single1, newpos1 )
            self.moveSingle( single2, newpos2 )


        return pair.single1, pair.single2


    def burstPair( self, pair ):

        single1, single2 = self.breakUpPair( pair )
        single1.initialize( self.t )
        single2.initialize( self.t )
        
        self.removeFromShellMatrix( pair )
        self.addToShellMatrix( single1 )
        self.addToShellMatrix( single2 )

        return single1, single2


    def formPairOrMulti( self, single, neighbors ):

        assert neighbors

        bursted = []

        # Try forming a Pair.
        if isinstance( neighbors[0], Single ):
            obj = self.formPair( single, neighbors[0], neighbors[1:] )
            if obj:
                return obj, neighbors[1:]


        # Then, a Multi.
        minShell = single.getMinSize() * ( 1.0 + self.MULTI_SHELL_FACTOR )
        # Good example why side-effect free programming is a good idea.  Then 
        # I would know for sure that formPair didn't do anything to neighbors. 
        neighborDists = self.objDistanceArray( single.pos, neighbors )
        neighbors = [ neighbors[i] for i in 
                      ( neighborDists <= minShell ).nonzero()[0] ]

        if not neighbors:
            return None, bursted

        closest = neighbors[0]

        if isinstance( closest, Single ):

            multi = self.createMulti()
            self.addToMulti( single, multi )
            self.removeFromShellMatrix( single )
            for neighbor in neighbors:
                self.addToMultiRecursive( neighbor, multi )

            multi.initialize( self.t )
            
            self.addToShellMatrix( multi )
            self.addMultiEvent( multi )

            return multi, bursted

        elif isinstance( closest, Multi ):

            multi = closest
            log.info( 'multi merge %s %s' % ( single, multi ) )

            self.removeFromShellMatrix( multi )

            self.addToMulti( single, multi )
            self.removeFromShellMatrix( single )
            for neighbor in neighbors[1:]:
                self.addToMultiRecursive( neighbor, multi )

            multi.initialize( self.t )

            self.addToShellMatrix( multi )
            self.updateEvent( self.t + multi.dt, multi )

            return multi, bursted


        assert False, 'do not reach here'


    def formPair( self, single1, single2, bursted ):

        #log.debug( 'trying to form %s' %
        #           'Pair( %s, %s )' % ( single1.particle, 
        #                                single2.particle ) )

        assert single1.isReset()
        assert single2.isReset()

        species1 = single1.particle.species
        species2 = single2.particle.species

        radius1 = species1.radius
        radius2 = species2.radius
        sigma = radius1 + radius2

        D1, D2 = species1.D, species2.D
        D12 = D1 + D2

        pairDistance = self.distance( single1.pos, single2.pos )
        r0 = pairDistance - sigma
        assert r0 >= 0, 'r0 (pair gap) between %s and %s = %g < 0' \
            % ( single1, single2, r0 )

        shellSize1 = pairDistance * D1 / D12 + radius1
        shellSize2 = pairDistance * D2 / D12 + radius2
        shellSizeMargin1 = radius1 * 2 #* self.SINGLE_SHELL_FACTOR
        shellSizeMargin2 = radius2 * 2 #* self.SINGLE_SHELL_FACTOR
        shellSizeWithMargin1 = shellSize1 + shellSizeMargin1
        shellSizeWithMargin2 = shellSize2 + shellSizeMargin2
        if shellSizeWithMargin1  >= shellSizeWithMargin2:
            minShellSize = shellSize1
            shellSizeMargin = shellSizeMargin1
        else:
            minShellSize = shellSize2
            shellSizeMargin = shellSizeMargin2

        # 1. Shell cannot be larger than max shell size or sim cell size.
        com = calculatePairCoM( single1.pos, single2.pos, D1, D2,
                                self.getWorldSize() )
        self.applyBoundary( com )
        minShellSizeWithMargin = minShellSize + shellSizeMargin
        maxShellSize = min( self.getMaxShellSize(),
                            r0 * 100 + sigma + shellSizeMargin )

        if minShellSizeWithMargin >= maxShellSize:
            log.debug( '%s not formed: minShellSize >= maxShellSize' %
                       ( 'Pair( %s, %s )' % ( single1.particle, 
                                              single2.particle ) ) )
            return None

        # Here, we have to take into account of the bursted Singles in
        # this step.  The simple check for closest below could miss
        # some of them, because sizes of these Singles for this
        # distance check has to include SINGLE_SHELL_FACTOR, while
        # these bursted objects have zero mobility radii.  This is not
        # beautiful, a cleaner framework may be possible.

        closest, closestShellDistance = DummySingle(), INF
        for b in bursted:
            if isinstance( b, Single ):
                d = self.distance( com, b.pos ) \
                    - b.getMinSize() * ( 1.0 + self.SINGLE_SHELL_FACTOR )
                if d < closestShellDistance:
                    closest, closestShellDistance = b, d

        if closestShellDistance <= minShellSizeWithMargin:
            log.debug( '%s not formed: squeezed by bursted neighbor %s' %
                       ( 'Pair( %s, %s )' % ( single1.particle, 
                                              single2.particle ), closest ) )
            return None


        c, d = self.getClosestObj( com, ignore=[ single1, single2 ] )
        if d < closestShellDistance:
            closest, closestShellDistance = c, d

        log.debug( 'Pair closest neighbor: %s %g, minShellWithMargin %g' %
                   ( closest, closestShellDistance, minShellSizeWithMargin ) )

        if isinstance( closest, Single ):

            D_closest = closest.particle.species.D
            D_tot = D_closest + D12
            closestDistance = self.distance( com, closest.pos )

            closestMinSize = closest.getMinSize()
            closestMinShell = closestMinSize * \
                ( self.SINGLE_SHELL_FACTOR + 1.0 )

            shellSize = min( ( D12 / D_tot ) *
                             ( closestDistance - minShellSize 
                               - closestMinSize ) + minShellSize,
                             closestDistance - closestMinShell,
                             closestShellDistance )

            shellSize /= SAFETY
            assert shellSize < closestShellDistance

        else:
            assert isinstance( closest, ( Pair, Multi, DummySingle ) )

            shellSize = closestShellDistance / SAFETY

        if shellSize <= minShellSizeWithMargin:
            log.debug( '%s not formed: squeezed by %s' %
                       ( 'Pair( %s, %s )' % ( single1.particle, 
                                              single2.particle ), closest ) )
            return None


        d1 = self.distance( com, single1.pos )
        d2 = self.distance( com, single2.pos )

        if shellSize < max( d1 + single1.getMinSize() *
                            ( 1.0 + self.SINGLE_SHELL_FACTOR ), \
                                d2 + single2.getMinSize() * \
                                ( 1.0 + self.SINGLE_SHELL_FACTOR ) ) * 1.3:
            log.debug( '%s not formed: singles are better' %
                       'Pair( %s, %s )' % ( single1.particle, 
                                            single2.particle ) )
            return None

        # 3. Ok, Pair makes sense.  Create one.
        shellSize = min( shellSize, maxShellSize )

        pair = self.createPair( single1, single2 )
        pair.setRadius( shellSize )

        self.removeFromShellMatrix( single1 )
        self.removeFromShellMatrix( single2 )
        self.addToShellMatrix( pair )

        pair.determineNextEvent( self.t )

        self.addPairEvent( pair )
        # single1 will be removed at the end of this step.
        self.removeEvent( single2 )

        assert closestShellDistance == INF or pair.radius < closestShellDistance
        assert pair.radius >= minShellSizeWithMargin
        assert pair.radius <= maxShellSize

        log.info( '%s, dt=%g, pairDistance=%g, shell=%g,' %
                  ( pair, pair.dt, pairDistance, pair.radius ) + 
                  'closest=%s, closestShellDistance=%g' %
                  ( closest, closestShellDistance ) )

        assert self.checkObj( pair )

        return pair
    


    def addToMultiRecursive( self, obj, multi ):
        
        if isinstance( obj, Single ):
            if obj.particle in multi.sim.particleList:  # Already in the Multi.
                return
            assert obj.isReset()
            
            self.addToMulti( obj, multi )
            self.removeFromShellMatrix( obj )
            self.removeEvent( obj )

            radius = obj.particle.species.radius *\
                ( 1.0 + self.MULTI_SHELL_FACTOR )
            neighbors = self.getNeighborsWithinRadiusNoSort( obj.pos, radius,
                                                             ignore=[obj,] )
            bursted = self.burstNonMultis( neighbors )
            neighborDists = self.objDistanceArray( obj.pos, bursted )
            neighbors = [ bursted[i] for i in 
                          ( neighborDists <= radius ).nonzero()[0] ]

            for obj in neighbors:
                self.addToMultiRecursive( obj, multi )

        elif isinstance( obj, Multi ):
            if not obj.sim.particleList[0] in multi.sim.particleList:
                self.mergeMultis( obj, multi )
                self.removeFromShellMatrix( obj )
                self.removeEvent( obj )
            else:
                log.debug( '%s already added. skipping.' % obj )
        else:
            assert False, 'do not reach here.'  # Pairs are bursted




    def addToMulti( self, single, multi ):
        log.info( 'adding %s to %s' % ( single, multi ) )

        shellSize = single.particle.species.radius * \
            ( 1.0 + self.MULTI_SHELL_FACTOR )
        multi.addParticle( single.particle )
        multi.addShell( single.pos, shellSize )


    '''
        merge multi1 into multi2
    '''
    def mergeMultis( self, multi1, multi2 ):

        log.info( 'merging %s to %s' % ( multi1, multi2 ) )

        assert not multi1.sim.particleList[0] in multi2.sim.particleList

        for i, particle in enumerate( multi1.sim.particleList ):
            
            # FIXME: shells should be renewed

            multi2.addParticle( particle )
            shell = multi1.shellList[i]
            multi2.addShell( shell.origin, shell.size )

        multi2.initialize( self.t )


    '''
    Find closest n shells.

    This method returns a tuple ( neighbors, distances ).
    '''

    # Sort yes/no.
    # Within radius yes/no.
    # With ignore list yes/no.

    # No radius (i.e. all), no ignore list.
    def getNeighborShells( self, pos, n=None ):
        neighbors, distances = self.sphereMatrix.getNeighbors( pos, n )

        if len( neighbors ) == 0:
            # Dummy
            return [( DummySingle(), 0 ),], [INF,]
        return neighbors, distances


    def getNeighborShellsNoSort( self, pos, n=None ):
        return self.sphereMatrix.getNeighborsNoSort( pos, n )


    # Within radius.
    def getNeighborShellsWithinRadius( self, pos, radius ):
        return self.sphereMatrix.getNeighborsWithinRadius( pos, radius )


    def getNeighborShellsWithinRadiusNoSort( self, pos, radius ):
        return self.sphereMatrix.getNeighborsWithinRadiusNoSort( pos, radius )


    # With radius and ignore list.
    def getNeighborsWithinRadius( self, pos, radius, ignore=[] ):
        shells, distances =\
            self.sphereMatrix.getNeighborsWithinRadius( pos, radius )

        neighbors = [ s[0] for s in shells if s[0] not in ignore ]
        neighbors = uniq( neighbors )
        return neighbors


    def getNeighborsWithinRadiusNoSort( self, pos, radius, ignore=[] ):
        shells, distances =\
            self.sphereMatrix.getNeighborsWithinRadiusNoSort( pos, radius )

        neighbors = uniq( [ s[0] for s in shells if s[0] not in ignore ] )
        return neighbors


    def getNeighbors( self, pos, radius=INF, ignore=[] ):
        shells, dists = self.sphereMatrix.getNeighbors( pos )

        seen = dict.fromkeys( ignore )
        neighbors = []
        distances = []

        for i, shell in enumerate( shells ):
            if not shell[0] in seen:
                seen[ shell ] = None
                neighbors.append(shell[0])
                distances.append(dists[i])
                if dists[i] > radius:
                    return neighbors, distances

        return neighbors + [DummySingle()], numpy.concatenate( [ distances,
                                                                 [INF] ] )


    def getClosestObj( self, pos, ignore=[] ):
        shells, distances = self.getNeighborShells( pos )

        for i, shell in enumerate( shells ):
            neighbor = shell[0]
            if neighbor not in ignore:
                return neighbor, distances[i]

        return DummySingle(), INF


    # Todo: do better for cylinders.
    def objDistance( self, pos, obj ):
        dists = numpy.zeros( len( obj.shellList ) )
        for i, shell in enumerate( obj.shellList ):
            dists[i] = self.distance( pos, shell.origin ) - shell.size
        return min( dists )

    """
    def Single.distanceTo( self, pos ):
        log.debug("Todo: better periodic boundary condition handling.")
        # Note: needs to work for multis as well.
        dists = numpy.zeros( len( self.shellList ) )
        for i, shell in enumerate( self.shellList ):
            dists[i] = shell.distanceTo( pos )
        return min( dists )
    """


    def objDistanceArray( self, pos, objs ):
        dists = numpy.array( [ self.objDistance( pos, obj ) for obj in objs ] )
        """
        log.debug("Todo: better periodic boundary condition handling.")
        dists = numpy.array( [ obj.distanceTo( pos ) for obj in objs ] )
        """
        return dists
            

    #
    # consistency checkers
    # Todo.
    #
    def checkObj( self, obj ):

        obj.check()

        allshells = [ ( obj, i ) for i in range( len( obj.shellList ) ) ]
        for i, shell in enumerate( obj.shellList ):

            closest, distance = self.getClosestObj( shell.origin,
                                                    ignore = [obj] )
            size = shell.size

            assert size <= self.getUserMaxShellSize(),\
                '%s shell size larger than user-set max shell size' % \
                str( ( obj, i ) )

            assert size <= self.getMaxShellSize(),\
                '%s shell size larger than simulator cell size / 2' % \
                str( ( obj, i ) )

            assert distance - size >= 0.0,\
                '%s overlaps with %s. (shell: %g, dist: %g, diff: %g.' \
                % ( str( obj ), str( closest ), size, distance,\
                        distance - size )

        return True


    def checkObjForAll( self ):

        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
            self.checkObj( obj )


    def checkEventStoichiometry( self ):

        population = 0
        for species in self.speciesList.values():
            population += species.pool.size

        eventPopulation = 0
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
            eventPopulation += obj.multiplicity

        if population != eventPopulation:
            raise RuntimeError, 'population %d != eventPopulation %d' %\
                  ( population, eventPopulation )


    def checkShellMatrix( self ):

        if self.worldSize != self.sphereMatrix.worldSize:
            raise RuntimeError,\
                'self.worldSize != self.sphereMatrix.worldSize'

        shellPopulation = 0
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
            shellPopulation += len( obj.shellList )

        if shellPopulation != self.sphereMatrix.size:
            raise RuntimeError,\
                'num shells != self.sphereMatrix.size'
        
        self.sphereMatrix.check()

        for k in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(k).getArg()
            for i in range( len( obj.shellList ) ):
                key = ( obj, i )
                pos, size = self.sphereMatrix.get( key )

                if ( obj.shellList[i].origin - pos ).any():
                    raise RuntimeError, \
                        '%s sphereMatrix positions consistency broken' % str( key )

                if obj.shellList[i].size != size:
                    raise RuntimeError, \
                        '%s sphereMatrix radii consistency broken' % str( key )



    def check( self ):

        ParticleSimulatorBase.check( self )

        assert self.scheduler.check()

        assert self.t >= 0.0
        assert self.dt >= 0.0

        self.checkShellMatrix()

        self.checkEventStoichiometry()
        
        self.checkObjForAll()




    #
    # methods for debugging.
    #

    def dumpScheduler( self ):
        scheduler = self.scheduler
        for i in range( scheduler.getSize() ):
            event = scheduler.getEventByIndex(i)
            print i, event.getTime(), event.getArg()

    def dump( self ):
        scheduler = self.scheduler
        for i in range( scheduler.getSize() ):
            event = scheduler.getEventByIndex(i)
            print i, event.getTime(), event.getArg(), event.getArg().pos


################ GRAVEYARD

    # Not used.
    """
    def getClosestShell( self, pos, ignore=[] ):
        neighbors, distances = self.getNeighborShells( pos )

        for i, neighbor in enumerate( neighbors ):
            if neighbor not in ignore:
                return neighbor, distances[i]

        return None, INF
    """

    # Not used.
    '''
    def getClosestNObjs( self, pos, n=1, ignore=[] ):

        neighbors, distances = self.getNeighborShells( pos, len( ignore ) + n )

        objs = []
        dists = []

        for i, neighbor in enumerate( neighbors ):
            if neighbor[0] not in ignore:
                objs += [neighbor[0]]
                dists += [distances[i]]
                if len( objs ) >= n:
                    return objs, dists

        return objs, dists
    '''

