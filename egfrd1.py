#!/usr/env python


import weakref

import math

import numpy
#import scipy
#import scipy.optimize


#from utils import *
from surface import *

from gfrdbase import *

from single import *
from pair import *
from multi import *


class Delegate( object ):

    def __init__( self, obj, method ):
        self.obj = weakref.proxy( obj )
        self.method = method

    def __call__( self, *arg ):
        return self.method( self.obj, *arg )


class EGFRDSimulator( ParticleSimulatorBase ):
    
    def __init__( self ):

        self.shellMatrix = ObjectMatrix()
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
        self.shellMatrix.setWorldSize( size )
        #self.sm2.setWorldSize( size )

    def setMatrixSize( self, size ):
        ParticleSimulatorBase.setMatrixSize( self, size )
        self.shellMatrix.setMatrixSize( size )
        #self.sm2.setMatrixSize( size )

    def getMatrixCellSize( self ):

        return self.shellMatrix.cellSize

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
        self.shellMatrix.clear()

        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                particle = Particle( species, index=i )
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





    def createSingle( self, particle ):

        rt = self.getReactionType1( particle.species )
        single = Single( particle, rt )
        single.initialize( self.t )

        return single


    def createPair( self, single1, single2 ):

        assert single1.dt == 0.0
        assert single2.dt == 0.0
        assert single1.getMobilityRadius() == 0.0
        assert single2.getMobilityRadius() == 0.0

        species1 = single1.particle.species
        species2 = single2.particle.species
        rt = self.reactionTypeMap2.get( ( species1, species2 ) )

        pair = Pair( single1, single2, rt, 
                     Delegate( self, EGFRDSimulator.distance ),
                     self.getWorldSize() )
        pair.initialize( self.t )

        return pair

    def createMulti( self ):

        multi = Multi( self )

        #multi.initialize( self.t )

        return multi


    def moveSingle( self, single, pos ):
        single.pos = pos
        self.updateOnParticleMatrix( single.particle, pos )


    def addToShellMatrix( self, obj ):
        for i, shell in enumerate( obj.shellList ):
            self.shellMatrix.add( ( obj, i ), shell.pos, shell.radius )


    def removeFromShellMatrix( self, obj ):
        for i in range( len( obj.shellList ) ):
            self.shellMatrix.remove( ( obj, i ) )


    def updateShellMatrix( self, obj ):
        for i, shell in enumerate( obj.shellList ):
            self.shellMatrix.update( ( obj, i ), shell.pos, shell.radius )


    def addEvent( self, t, func, arg ):

        return self.scheduler.addEvent( t, func, arg )

    def addSingleEvent( self, single ):
        eventID = self.addEvent( self.t + single.dt, 
                                 Delegate( self, EGFRDSimulator.fireSingle ), 
                                 single )
        single.eventID = eventID

    def addPairEvent( self, pair ):

        eventID = self.addEvent( self.t + pair.dt, 
                                 Delegate( self, EGFRDSimulator.firePair ), 
                                 pair )
        pair.eventID = eventID

    def addMultiEvent( self, multi ):

        eventID = self.addEvent( self.t + multi.dt, 
                                 Delegate( self, EGFRDSimulator.fireMulti ), 
                                 multi )
        multi.eventID = eventID


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



    def fireSingleReaction( self, single ):

        reactantSpecies = single.particle.species
        oldpos = single.particle.pos.copy()
        
        rt = single.drawReactionType()

        if len( rt.products ) == 0:
            
            self.removeParticle( single.particle )

            self.lastReaction = Reaction( rt, [single.particle], [] )

            
        elif len( rt.products ) == 1:
            
            productSpecies = rt.products[0]

            if reactantSpecies.radius < productSpecies.radius:
                self.clearVolume( oldpos, productSpecies.radius )

            if not self.checkOverlap( oldpos, productSpecies.radius,
                                      ignore = [ single.particle, ] ):
                log.info( 'no space for product particle.' )
                raise NoSpace()

            self.removeParticle( single.particle )
            newparticle = self.createParticle( productSpecies, oldpos )
            newsingle = self.createSingle( newparticle )
            self.addToShellMatrix( newsingle )
            self.addSingleEvent( newsingle )

            self.lastReaction = Reaction( rt, [single.particle], [newparticle] )

            log.info( 'product; %s' % str( newsingle ) )

            
        elif len( rt.products ) == 2:
            
            productSpecies1 = rt.products[0]
            productSpecies2 = rt.products[1]
            
            D1 = productSpecies1.D
            D2 = productSpecies2.D
            D12 = D1 + D2
            
            particleRadius1 = productSpecies1.radius
            particleRadius2 = productSpecies2.radius
            particleRadius12 = particleRadius1 + particleRadius2

            # clean up space.
            rad = max( particleRadius12 * ( D1 / D12 ) + particleRadius1,
                       particleRadius12 * ( D2 / D12 ) + particleRadius2 )

            self.clearVolume( oldpos, rad )

            for _ in range( 100 ):
                unitVector = randomUnitVector()
                vector = unitVector * particleRadius12 * ( 1.0 + 1e-7 )
            
                # place particles according to the ratio D1:D2
                # this way, species with D=0 doesn't move.
                # FIXME: what if D1 == D2 == 0?

                while 1:
                    newpos1 = oldpos + vector * ( D1 / D12 )
                    newpos2 = oldpos - vector * ( D2 / D12 )
                    self.applyBoundary( newpos1 )
                    self.applyBoundary( newpos2 )

                    if self.distance( newpos1, newpos2 ) >= particleRadius12:
                        break

                    vector *= 1.0 + 1e-7


                # accept the new positions if there is enough space.
                if ( self.checkOverlap( newpos1, particleRadius1,
                                        ignore = [ single.particle, ] ) and
                     self.checkOverlap( newpos2, particleRadius2,
                                        ignore = [ single.particle, ] ) ):
                    break
            else:
                log.info( 'no space for product particles.' )
                raise NoSpace()

            self.removeParticle( single.particle )

            particle1 = self.createParticle( productSpecies1, newpos1 )
            particle2 = self.createParticle( productSpecies2, newpos2 )
            newsingle1 = self.createSingle( particle1 )
            newsingle2 = self.createSingle( particle2 )

            self.addToShellMatrix( newsingle1 )
            self.addToShellMatrix( newsingle2 )
            self.addSingleEvent( newsingle1 )
            self.addSingleEvent( newsingle2 )

            self.lastReaction = Reaction( rt, [single.particle], 
                                          [particle1, particle2] )

            log.info( 'products; %s %s' % 
                      ( str( newsingle1 ), str( newsingle2 ) ) )

        else:
            raise RuntimeError, 'num products >= 3 not supported.'

        self.reactionEvents += 1


    def propagateSingle( self, single, r ):

        assert abs( single.dt + single.lastTime - self.t ) <= 1e-18 * self.t
        
        displacement = randomVector( r )

        assert abs( length( displacement ) - r ) <= 1e-15 * r
            
        newpos = single.particle.pos + displacement
        self.applyBoundary( newpos )
            
        assert self.checkOverlap( newpos, single.getMinRadius(),
                                  ignore = [ single.particle, ] )

        self.moveSingle( single, newpos )

        single.initialize( self.t )

        self.updateShellMatrix( single )



    def fireSingle( self, single ):

        # Reaction.
        if single.eventType == EventType.REACTION:

            log.info( 'single reaction %s' % str( single ) )
            r = single.drawR( single.dt )

            self.propagateSingle( single, r )

            try:
                self.removeFromShellMatrix( single )
                self.fireSingleReaction( single )
            except NoSpace:
                log.info( 'single reaction; placing product failed.' )
                self.addToShellMatrix( single )
                self.rejectedMoves += 1
                single.reset()
                return single.dt

            single.dt = -INF  # remove this Single from the Scheduler
            return single.dt

        # Propagate, if not reaction.

        # Handle immobile case first.
        if single.getD() == 0:
            # no propagation, just calculate next reaction time.
            single.determineNextEvent( self.t ) 
            return single.dt
        
        # Propagate this particle to the exit point on the shell.
        
        self.propagateSingle( single, single.getMobilityRadius() )

        # (2) Clear volume.

        minShell = single.getMinRadius() * ( 1.0 + self.SINGLE_SHELL_FACTOR )

        closeNeighbors, distances = self.getNeighbors( single.pos, minShell,
                                                       ignore=[single,] )

        # This is a bit tricky, but the last one in closeNeighbors
        # is the closest object to this Single.
        # getNeighbors() returns closeNeighbors within minShell *plus* one.
        closest = closeNeighbors.pop()
        closestShellDistance = distances[-1]

        bursted = []
        
        if closeNeighbors:
            bursted = self.burstNonMultis( closeNeighbors )
            obj, b = self.formPairOrMulti( single, bursted )
            bursted.extend( b )

            if obj:
                single.dt = -INF # remove by rescheduling to past.
                return single.dt

            # if nothing was formed, recheck closest and restore shells.
            closest, closestShellDistance = \
                self.getClosestObj( single.pos, ignore = [ single, ] )

        self.updateSingle( single, closest, closestShellDistance )

        bursted = uniq( bursted )
        burstedSingles = [ s for s in bursted if isinstance( s, Single ) ]
        self.restoreSingleShells( burstedSingles )
            
        log.info( 'single shell %g dt %g.' % ( single.radius, single.dt ) )

        return single.dt


    def restoreSingleShells( self, singles ):

        for single in singles:
            assert single.isReset()
            c, d = self.getClosestObj( single.pos, ignore = [single,] )

            self.updateSingle( single, c, d )
            self.updateEvent( self.t + single.dt, single )
            log.debug( 'restore shell %s %g dt %g closest %s %g' %
                       ( single, single.radius, single.dt, c, d ) )


    def calculateSingleShellSize( self, single, closest, 
                                  distance, shellDistance ):

        minRadius1 = single.getMinRadius()
        D1 = single.getD()

        if D1 == 0:
            return minRadius1

        D2 = closest.getD()
        minRadius2 = closest.getMinRadius()
        minRadius12 = minRadius1 + minRadius2
        sqrtD1 = math.sqrt( D1 )
            
        shellSize = min( sqrtD1 / ( sqrtD1 + math.sqrt( D2 ) )
                         * ( distance - minRadius12 ) + minRadius1,
                         shellDistance )
        shellSize /= SAFETY
        shellSize = max( shellSize, minRadius1 ) # not smaller than the radius

        return shellSize


    def updateSingle( self, single, closest, distanceToShell ): 

        if isinstance( closest, Single ):
            distanceToClosest = self.distance( single.pos, closest.pos )
            shellSize = self.calculateSingleShellSize( single, closest, 
                                                       distanceToClosest,
                                                       distanceToShell )
        else:  # Pair or Multi
            shellSize = distanceToShell / SAFETY
            shellSize = max( shellSize, single.getMinRadius() )

        shellSize = min( shellSize, self.getMaxShellSize() )

        single.setRadius( shellSize )
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
                self.applyBoundary( newCoM )

                self.removeParticle( particle1 )
                self.removeParticle( particle2 )

                particle = self.createParticle( species3, newCoM )
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


    def burstSingle( self, single ):

        assert self.t >= single.lastTime
        assert self.t <= single.lastTime + single.dt
        assert single.radius >= single.getMinRadius()

        dt = self.t - single.lastTime

        particleRadius = single.particle.species.radius
        oldpos = single.particle.pos # .copy()

        r = single.drawR( dt )
        displacement = randomVector( r )

        newpos = oldpos + displacement

        self.applyBoundary( newpos )

        assert self.distance( newpos, oldpos ) <= single.getMobilityRadius()
        assert self.distance( newpos, oldpos ) - r <= r * 1e-6
        assert self.checkOverlap( newpos, particleRadius,\
                                  ignore = [ single.particle, ] )

        self.moveSingle( single, newpos )

        single.initialize( self.t )
        self.updateShellMatrix( single )
        self.updateEvent( self.t, single )


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
        minShell = single.getMinRadius() * ( 1.0 + self.MULTI_SHELL_FACTOR )
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
                    - b.getMinRadius() * ( 1.0 + self.SINGLE_SHELL_FACTOR )
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

            closestMinRadius = closest.getMinRadius()
            closestMinShell = closestMinRadius * \
                ( self.SINGLE_SHELL_FACTOR + 1.0 )

            shellSize = min( ( D12 / D_tot ) *
                             ( closestDistance - minShellSize 
                               - closestMinRadius ) + minShellSize,
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

        if shellSize < max( d1 + single1.getMinRadius() *
                            ( 1.0 + self.SINGLE_SHELL_FACTOR ), \
                                d2 + single2.getMinRadius() * \
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
            multi2.addShell( shell.pos, shell.radius )

        multi2.initialize( self.t )


    '''
    Find closest n shells.

    This method returns a tuple ( neighbors, distances ).
    '''

    def getNeighborShells( self, pos, n=None ):

        neighbors, distances = self.shellMatrix.getNeighbors( pos, n )

        if len( neighbors ) == 0:
            return [( DummySingle(), 0 ),], [INF,]
        return neighbors, distances


    def getNeighborShellsNoSort( self, pos, n=None ):

        return self.shellMatrix.getNeighborsNoSort( pos, n )


    def getNeighborShellsWithinRadius( self, pos, radius ):
        return self.shellMatrix.getNeighborsWithinRadius( pos, radius )


    def getNeighborShellsWithinRadiusNoSort( self, pos, radius ):
        return self.shellMatrix.getNeighborsWithinRadiusNoSort( pos, radius )


    def getNeighborsWithinRadius( self, pos, radius, ignore=[] ):

        shells, distances =\
            self.shellMatrix.getNeighborsWithinRadius( pos, radius )

        neighbors = [ s[0] for s in shells if s[0] not in ignore ]
        neighbors = uniq( neighbors )

        return neighbors


    def getNeighborsWithinRadiusNoSort( self, pos, radius, ignore=[] ):

        shells, distances =\
            self.shellMatrix.getNeighborsWithinRadiusNoSort( pos, radius )

        neighbors = uniq( [ s[0] for s in shells if s[0] not in ignore ] )

        return neighbors


    def getNeighbors( self, pos, radius=INF, ignore=[] ):

        shells, dists = self.shellMatrix.getNeighbors( pos )

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

    def getClosestShell( self, pos, ignore=[] ):

        neighbors, distances = self.getNeighborShells( pos )

        for i, neighbor in enumerate( neighbors ):
            if neighbor not in ignore:
                return neighbor, distances[i]

        return None, INF


    def getClosestObj( self, pos, ignore=[] ):

        shells, distances = self.getNeighborShells( pos )

        for i, shell in enumerate( shells ):
            neighbor = shell[0]
            if neighbor not in ignore:
                return neighbor, distances[i]

        return DummySingle(), INF


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

    def objDistance( self, pos, obj ):
        
        dists = numpy.zeros( len( obj.shellList ) )
        for i, shell in enumerate( obj.shellList ):
            dists[i] = self.distance( pos, shell.pos ) - shell.radius

        return min( dists )

    def objDistanceArray( self, pos, objs ):

        dists = numpy.array( [ self.objDistance( pos, obj ) for obj in objs ] )
        return dists
            

    #
    # consistency checkers
    #

    
    def checkObj( self, obj ):

        obj.check()

        allshells = [ ( obj, i ) for i in range( len( obj.shellList ) ) ]
        for i, shell in enumerate( obj.shellList ):

            closest, distance = self.getClosestObj( shell.pos,
                                                    ignore = [obj] )
            radius = shell.radius

            assert radius <= self.getUserMaxShellSize(),\
                '%s shell size larger than user-set max shell size' % \
                str( ( obj, i ) )

            assert radius <= self.getMaxShellSize(),\
                '%s shell size larger than simulator cell size / 2' % \
                str( ( obj, i ) )

            assert distance - radius >= 0.0,\
                '%s overlaps with %s. (shell: %g, dist: %g, diff: %g.' \
                % ( str( obj ), str( closest ), radius, distance,\
                        distance - radius )

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

        if self.worldSize != self.shellMatrix.worldSize:
            raise RuntimeError,\
                'self.worldSize != self.shellMatrix.worldSize'

        shellPopulation = 0
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
            shellPopulation += len( obj.shellList )

        if shellPopulation != self.shellMatrix.size:
            raise RuntimeError,\
                'num shells != self.shellMatrix.size'
        
        self.shellMatrix.check()

        for k in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(k).getArg()
            for i in range( len( obj.shellList ) ):
                key = ( obj, i )
                pos, radius = self.shellMatrix.get( key )

                if ( obj.shellList[i].pos - pos ).any():
                    raise RuntimeError, \
                        '%s shellMatrix positions consistency broken' % str( key )

                if obj.shellList[i].radius != radius:
                    raise RuntimeError, \
                        '%s shellMatrix radii consistency broken' % str( key )



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



