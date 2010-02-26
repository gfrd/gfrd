#!/usr/env python
import math
import numpy

from gfrdbase import *

from single import *
from pair import *
from multi import *
from shape import *

from log import *
from utils import *


class EGFRDSimulator( ParticleSimulatorBase ):
    """
    Return a eGFRD simulator which contains the geometry and particles and 
    carries out the simulation.

    A cubic simulation unit cell from [0,0,0] to [L,L,L] is created, if L is 
    the argument.
    """

    def __init__( self, worldSize ):
        self.sphereMatrix = SphereMatrix()
        self.cylinderMatrix = CylinderMatrix()
        self.objectMatrices = [ self.sphereMatrix, self.cylinderMatrix ]
        #self.sm2 = pObjectMatrix()

        ParticleSimulatorBase.__init__( self, worldSize )

        self.isDirty = True
        self.scheduler = EventScheduler()

        #self.smallT = 1e-8  # FIXME: is this ok?
        self.userMaxShellSize = INF
        self.reset()


    def setWorldSize( self, size ):
        """Called from ParticleSimulatorBase.__init__()

        """
        if isinstance( size, list ) or isinstance( size, tuple ):
            size = numpy.array( size )

        ParticleSimulatorBase.setWorldSize( self, size )

        self.sphereMatrix.setWorldSize( size )
        self.cylinderMatrix.setWorldSize( size )


    def setUserMaxShellSize( self, size ):
        self.userMaxShellSize = size


    def getUserMaxShellSize( self ):
        return self.userMaxShellSize


    def getMaxShellSize( self ):
        return min( self.getMatrixCellSize() * .5 / SAFETY,
                    self.userMaxShellSize )


    def getNextTime( self ):
        if self.scheduler.getSize() == 0:
            return self.t

        return self.scheduler.getTopTime()


    def reset( self ):
        self.t = 0.0
        self.dt = 0.0
        self.stepCounter = 0
        self.zeroSteps = 0
        self.rejectedMoves = 0
        self.reactionEvents = 0
        self.unbindingEvents = 0
        self.bindingEvents = 0
        self.lastEvent = None
        self.lastReaction = None
        self.isDirty = True
        #self.initialize()

        self.single_steps = {EventType.ESCAPE:0, EventType.REACTION:0}
        self.pair_steps = {0:0,1:0,2:0,3:0}
        self.multi_steps = {EventType.ESCAPE:0, EventType.REACTION:0, 2:0}

    def initialize( self, a=None ):
        ParticleSimulatorBase.initialize( self )
        self.setAllRepulsive()
        self.scheduler.clear()
        self.sphereMatrix.clear()
        self.cylinderMatrix.clear()

        # Get particle data from particleMatrix, because we need to know the 
        # surface they are on.
        particles = self.particleMatrix.getAll( )
        for particle in particles:
            single = self.createSingle( particle )

        nParticles = sum( [ s.pool.size for s in self.speciesList.values() ] )
        assert nParticles == len( particles )
        self.isDirty = False


    def step( self ):
        if self.isDirty:
            self.initialize( 1 )
        self.lastReaction = None
            
        '''
        if self.stepCounter % 100 == 0:
            self.check()
        '''
        
        self.stepCounter += 1
        event = self.scheduler.getTopEvent()
        self.t, self.lastEvent = event.getTime(), event.getArg()

        log.info( '\n\n    %d: t = %s dt = %s\n    %s.\n    '
                  'Reactions = %d, unbindings = %d, bindings = %d, '
                  'rejectedmoves = %d.' %
                  ( self.stepCounter, repr( self.t ), repr( self.dt ), 
                    self.lastEvent, self.reactionEvents, self.unbindingEvents,
                    self.bindingEvents, self.rejectedMoves ) )

        self.scheduler.step()

        if self.particleMatrix.size == 0:
            raise RuntimeError( '0 particles left' )

        nextTime = self.scheduler.getTopTime()
        self.dt = nextTime - self.t

        # Assert if not too many successive dt == 0 steps occur.
        if self.dt == 0.0:
            self.zeroSteps += 1
            if self.zeroSteps >= max( self.scheduler.getSize() * 3, 10 ):
                raise RuntimeError( 'Too many dt = zero steps. '
                                    'Simulator halted?' )
        else:
            self.zeroSteps = 0
        assert self.scheduler.getSize() != 0


    def stop( self, t ):
        '''Request for stop at time t.

        '''
        log.info( '    Stop at %s.' % repr( t ) )

        nextEventTime = self.scheduler.getTopEvent().getTime()
        if self.t == t:
            # Todo. Why return here?
            return
        if t >= nextEventTime:
            # There is an event in the scheduler that needs to be handled 
            # first.
            raise RuntimeError( 'Stop time %s >= next event time %s. ' % 
                                ( repr( t ), repr( nextEventTime ) ) )
        if t < self.t:
            raise RuntimeError( 'Stop time %s < current time %s.' %
                                ( repr( t ), repr( self.t ) ) )

        self.t = t
        scheduler = self.scheduler
        nonSingleList = []

        # First burst all Singles.
        for i in range( scheduler.getSize() ):
            obj = scheduler.getEventByIndex( i ).getArg()
            if isinstance( obj, Pair ) or isinstance( obj, Multi ):
                nonSingleList.append( obj )
            elif isinstance( obj, Single ):
                log.debug( '        *Burst %s, lastTime = %.3g' %
                           ( str( obj ), obj.lastTime ) )
                self.burstSingle( obj )
            else:
                assert False, 'do not reach here'


        # then burst all Pairs and Multis.
        log.debug( '        *Burst %s' % nonSingleList )
        self.burstObjs( nonSingleList )
        self.dt = 0.0

    
    ##########################################################################
    # Part of the burst methods.

    def clearVolume( self, pos, radius, ignore=[] ):
        neighbors = self.getNeighborsWithinRadiusNoSort( pos, radius, ignore )
        return self.burstObjs( neighbors )


    def burstNonMultis( self, neighbors ):
        bursted = []
        for obj in neighbors:
            if not ( isinstance( obj, Multi ) or isinstance( obj, Surface ) ):
                b = self.burstObj( obj )
                bursted.extend( b )
            else:
                bursted.append( obj )
        return bursted


    def burstObjs( self, objs ):
        bursted = []
        for obj in objs:
            b = self.burstObj( obj )
            bursted.extend( b )

        return bursted


    def burstObj( self, obj ):
        log.info( '    *Bursting\n'
                  '        %s.' % ( str( obj ) ) )
        if isinstance( obj, Single ):
            obj = self.burstSingle( obj )
            return [ obj, ]
        elif isinstance( obj, Pair ):  # Pair
            single1, single2 = self.burstPair( obj )
            self.removeEvent( obj )
            self.addSingleEvent( single1 )
            self.addSingleEvent( single2 )
            return [ single1, single2 ]
        elif isinstance( obj, Multi ): # Multi
            bursted = self.burstMulti( obj )
            self.removeEvent( obj )
            return bursted
        else:
            raise NotImplementedError()


    ##########################################################################
    # Event scheduler stuff.

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


    def addEvent( self, t, func, arg ):
        return self.scheduler.addEvent( t, func, arg )


    def removeEvent( self, event ):
        self.scheduler.removeEvent( event.eventID )


    def updateEvent( self, t, event ):
        self.scheduler.updateEventTime( event.eventID, t )


    ##########################################################################
    # Shell matrix stuff.

    def setMatrixSize( self, size ):
        ParticleSimulatorBase.setMatrixSize( self, size )
        self.sphereMatrix.setMatrixSize( size )
        self.cylinderMatrix.setMatrixSize( size )


    def getMatrixCellSize( self ):
        assert self.sphereMatrix.cellSize == self.cylinderMatrix.cellSize
        return self.sphereMatrix.cellSize


    def addToShellMatrix( self, obj ):
        """Add all shells from obj.shellList to one of the shell matrices.

        """
        for i, shell in enumerate( obj.shellList ):
            key = ( obj, i )
            if isinstance( shell, Sphere ):
                self.sphereMatrix.add( key, shell.origin, shell.radius )
            elif isinstance( shell, Cylinder ):
                self.cylinderMatrix.add( key, shell )
            else: raise KeyError( 'Objecttype does not exit' )


    def updateShellMatrix( self, obj ):
        for i, shell in enumerate( obj.shellList ):
            key = ( obj, i )
            if isinstance( shell, Sphere ):
                self.sphereMatrix.update( key, shell.origin, shell.radius )
            elif isinstance( shell, Cylinder ):
                self.cylinderMatrix.update( key, shell )
            else: raise KeyError( 'Objecttype does not exit' )


    def removeFromShellMatrix( self, obj ):
        for i, shell in enumerate( obj.shellList ):
            key = ( obj, i )
            if isinstance( shell, Sphere ):
                self.sphereMatrix.remove( key )
            elif isinstance( shell, Cylinder ):
                self.cylinderMatrix.remove( key )
            else: raise KeyError( 'Objecttype does not exit' )


    ##########################################################################
    # Find nearest neighbors.
    #
    # Most of these methods (also see cObjectmatrix.py) return a tuple 
    # (neighbors, distances).
    #
    # Options:
    #   * Sort yes/no.
    #   * Within radius yes/no.
    #   * With ignore list yes/no.

    def getClosestObj( self, pos, ignore=[] ):
        """getClosestObject.


        No radius (i.e. all), no ignore list, keys not extracted.

        Don't be smart and do:
            keys, distances, _, _ = self.getNeighbors( pos, ignore )
            return keys[0], distance[0]
        since that would unneccesarily iterate one more time over all 
        objects in self.getNeighbors. 

        """
        closestObject = DummySingle()
        closestDistance = INF

        for objectMatrix in self.objectMatrices:
            keys, distances = objectMatrix.getNeighbors( pos )

            # Don't be clever and think you can do:
            #    keys, distances = objectMatrix.getNeighbors( pos, 1 )
            # Because then you are ignoring the ignore list.

            for i, key in enumerate( keys ):
                object = key[0]
                if object not in ignore and distances[i] < closestDistance:
                    closestObject, closestDistance = object, distances[i]
                    # Found yet a closer single. Break out of inner for loop 
                    # and check other objectMatrices.
                    break   

        # Surface detection.
        distanceToSurface, closestSurface = self.getClosestSurface( pos, 
                                                                    ignore ) 
        if distanceToSurface < closestDistance:
            closestDistance = distanceToSurface
            closestObject = closestSurface

        return closestObject, closestDistance


    def getNeighborsWithinRadiusNoSort( self, pos, radius, ignore=[] ):
        """Only returns neighbors, not the distances towards their shells. Can 
        for example be used to try to clear all objects from a certain volume.

        """
        neighbors = []
        for objectMatrix in self.objectMatrices:
            keys, _ = \
                objectMatrix.getNeighborsWithinRadiusNoSort( pos, radius )
            # object = key[0]
            # Since an object can have more than 1 shell (multis for example), 
            # and for each shell there is an entry in the objectMatrix, we 
            # make sure each object occurs only once in the returned list 
            # here.
            neighbors.extend( uniq( [ key[0] for key in keys 
                                      if key[0] not in ignore ] ) )
        return neighbors


    def getNeighbors( self, pos, radius=INF, ignore=[] ):
        """Return a sorted list of all objects and the distances towards their 
        *shell* if that shell is within 'radius' of 'pos', plus the closest 
        object and the distance towards its shell that is just outside of 
        radius. 

        Read that sentence again, it's all in there.

        Only called once from fireSingle.

        """
        # Diccionaries are more efficient for lookups.
        seen = dict.fromkeys( ignore )
        neighbors = []
        distances = []

        closestObject = DummySingle()
        closestDistance = INF

        for objectMatrix in self.objectMatrices:
            # Don't name this distances, it's already defined.
            keys, dists = objectMatrix.getNeighbors( pos )

            for i, key in enumerate( keys ):
                distance = dists[i]
                object = key[0]
                if object.dt == 0:
                    # Make bursted singles look bigger. Right now the size of 
                    # their shell is only particle.radius (not yet multiplied 
                    # by MINIMAL_SINGLE_RADIUS_FACTOR), and that might cause 
                    # the simulation might to come to a halt. Single 1 bursts 
                    # single 2 and forms a shell, next step single 2 bursts 
                    # single 1 and forms a shell, etc.
                    distance -= object.radius * \
                                ( MINIMAL_SINGLE_RADIUS_FACTOR - 1 )

                if not object in seen:
                    # Since an object can have more than 1 shell (multis for 
                    # example), and for each shell there is an entry in the 
                    # objectMatrix, we signal here that we have already seen 
                    # and added this object. Since the 
                    # objectMatrix.getNeighbors returns a sorted list, this 
                    # works.
                    #
                    # Todo. This used to be a bug, right?: seen[ key ] = None
                    #
                    # Any value other then 'None' would be fine here also.
                    seen[ object ] = None
                    if distance > radius:
                        # This is an object (the first for this objectMatrix) 
                        # that has a shell that is more than radius away from 
                        # pos. If it is closer than the closest such one we 
                        # found so far: store it. Always break out of the 
                        # inner for loop now and check the other 
                        # objectMatrices.
                        if  distance < closestDistance:
                            closestObject = object
                            closestDistance = distance
                            break
                        else:
                            break # Just to be sure you get the point.
                    else:
                        # This is an object that has a shell that is within 
                        # 'radius' of pos, so add it to the neighbors list.
                        neighbors.append( object )
                        distances.append( distance )

        # Surface detection.
        distanceToSurface, closestSurface = self.getClosestSurface( pos, 
                                                                    ignore ) 
        if distanceToSurface < closestDistance:
            closestDistance = distanceToSurface
            closestObject = closestSurface

        return neighbors, distances, closestObject, closestDistance


    def objDistance( self, pos, obj ):
        """objDistance and objDistanceArray can be used to get information on 
        the distance from a certain point to some *previously found* objects. 
        So no call to the objectmatrices is made here.

        Possible improvement: get this info from objectMatrix using 
        objectMatrix.get().

        """
        if isinstance( obj, Surface ):
            posTransposed = cyclicTranspose( pos, obj.origin, self.worldSize )
            return obj.signedDistanceTo( posTransposed )
        else:
            dists = numpy.zeros( len( obj.shellList ) )
            # A multi (and possibly singles and pairs as well) has several 
            # shells, loop over them.
            for i, shell in enumerate( obj.shellList ):
                posTransposed = cyclicTranspose( pos, shell.origin,
                                                 self.worldSize )
                dists[i] = shell.signedDistanceTo( posTransposed )
            return min( dists )


    def objDistanceArray( self, pos, objs ):
        dists = numpy.array( [ self.objDistance( pos, obj ) for obj in objs ] )
        return dists
            

    ##########################################################################
    # consistency checkers
    # Todo.

    def checkObj( self, obj ):
        #log.debug( '        *checkObj: %s.' % ( obj ) )
        obj.check()

        allshells = [ ( obj, i ) for i in range( len( obj.shellList ) ) ]
        for i, shell in enumerate( obj.shellList ):
            closest, distance = self.getClosestObj( shell.origin,
                                                    ignore=[ obj, ] )
            radius = shell.radius


            assert radius <= self.getUserMaxShellSize(), \
                   '%s shell radius larger than user-set max shell radius' % \
                   str( ( obj, i ) )
            assert radius <= self.getMaxShellSize(), \
                   '%s shell radius larger than simulator cell radius / 2. ' \
                   '%.3g vs %.3g' % \
                   ( str( ( obj, i ) ), radius, self.getMaxShellSize() )

            # Todo.
            break
            
            if not ( isinstance( obj, InteractionSingle ) and \
                     isinstance( closest, Surface ) ):
                assert distance - radius >= 0.0, \
                       '%s overlaps with %s. shell: %.3g, dist: %.3g, ' \
                       'diff: %.3g.' % ( str( obj ), str( closest ), radius,
                                         distance, distance - radius )
        return True


    def checkObjForAll( self ):
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex( i ).getArg()
            self.checkObj( obj )


    def checkEventStoichiometry( self ):
        population = 0
        for species in self.speciesList.values():
            population += species.pool.size

        eventPopulation = 0
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex( i ).getArg()
            eventPopulation += obj.multiplicity

        if population != eventPopulation:
            raise RuntimeError( 'population %d != eventPopulation %d' %
                                ( population, eventPopulation ) )


    def checkShellMatrix( self ):
        if self.worldSize != self.sphereMatrix.worldSize:
            raise RuntimeError( 'self.worldSize != '
                                'self.sphereMatrix.worldSize' )

        shellPopulation = 0
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex( i ).getArg()
            shellPopulation += len( obj.shellList )
        
        matrixPopulation = sum( matrix.size for matrix in self.objectMatrices )
        if shellPopulation != matrixPopulation:
            raise RuntimeError( 'num shells %.3g != matrixPopulation %.3g' % 
                                ( shellPopulation, matrixPopulation ) )
        
        self.sphereMatrix.check()

        for k in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex( k ).getArg()
            for i in range( len( obj.shellList ) ):
                key = ( obj, i )
                if isinstance( key[0].shellList[0], Cylinder ):
                    # Todo.
                    break

                pos, radius = self.sphereMatrix.get( key )

                if ( obj.shellList[i].origin - pos ).any():
                    raise RuntimeError( '%s sphereMatrix positions consistency ' 
                                        'broken' % str( key ) )

                if obj.shellList[i].radius != radius:
                    raise RuntimeError( '%s sphereMatrix radii consistency '
                                        'broken' % str( key ) )


    def check( self ):
        ParticleSimulatorBase.check( self )

        assert self.scheduler.check()

        assert self.t >= 0.0
        assert self.dt >= 0.0

        self.checkShellMatrix()

        self.checkEventStoichiometry()
        
        self.checkObjForAll()


    ##########################################################################
    # Methods for debugging.

    def dumpScheduler( self ):
        scheduler = self.scheduler
        for i in range( scheduler.getSize() ):
            event = scheduler.getEventByIndex( i )
            print i, event.getTime(), event.getArg()


    def dump( self ):
        scheduler = self.scheduler
        for i in range( scheduler.getSize() ):
            event = scheduler.getEventByIndex( i )
            print i, event.getTime(), event.getArg(), event.getArg().pos


    ##########################################################################
    # Methods for Singles.

    def createSingle( self, particle ):
        # Poissonian reactions.
        reactionType = self.getReactionType1( particle.species )
        # The type of defaultSingle depends on the surface this particle is 
        # on. Either SphericalSingle, PlanarSurfaceSingle, or 
        # CylindricalSurfaceSingle.
        single = particle.surface.defaultSingle( particle, reactionType )
        single.initialize( self.t )
        self.addToShellMatrix( single )
        self.addSingleEvent( single )
        return single


    def fireSingle( self, single ):
        """Exectute unimolecular reaction, or escape event for single. Return 
        new dt to the scheduler to reschedule the event, or return -INF to 
        remove it. 

        """
        assert abs( single.dt + single.lastTime - self.t ) <= 1e-18 * self.t

        # 0. Reaction.
        if single.eventType == EventType.REACTION:
            self.propagateSingle( single )
            self.single_steps[single.eventType] += 1
            log.info( '    ' + str( single.eventType ) )

            try:
                self.fireSingleReaction( single )
            except NoSpace:
                log.info( '    single reaction; placing product failed.' )
                self.addToShellMatrix( single )
                self.rejectedMoves += 1
                single.reset()
                return single.dt

            # remove old single event from the Scheduler
            single.dt = -INF
            return single.dt

        # 1. Escape or Interaction.
        # Handle immobile case first.
        if single.getD() == 0:
            # No propagation, just calculate next reaction time.
            # (!) Not a pure function.
            single.dt, single.eventType, single.activeDomain = \
                single.determineNextEvent( ) 
            single.lastTime = self.t
            return single.dt

        # Propagate if dt != 0.
        if single.dt != 0.0:
            # All Singles, including NonInteractionSingles, should call 
            # drawEventType, even if it is already known that it is going to 
            # be an escape. This is because drawEventType also sets the escape 
            # flag in the active domain, which is used when calling 
            # propagateSingle and signals that in this domain the exitpoint is 
            # already known. Read that again.
            single.eventType = single.activeDomain.drawEventType( single.dt )

            # Propagate this particle to the exit point on the shell.
            self.propagateSingle( single )

        self.single_steps[single.eventType] += 1

        if single.eventType == EventType.REACTION:
            # Reactions and interactions make use of the same eventType, but 
            # we know this must be an interaction because reactions have been 
            # taken care of above.
            log.info( '    INTERACTION' )
            try:
                self.fireSingleReaction( single, True )
            except NoSpace:
                log.info( '    Single interaction; placing product failed.' )
                self.addToShellMatrix( single )
                self.rejectedMoves += 1
                single.reset()
                return single.dt
        else:
            log.info( '    ESCAPE' )


        if isinstance( single, InteractionSingle ):
            # A new single is created, either in propagateSingle or in 
            # fireSingleReaction. Remove the old single from scheduler. Both 
            # if eventType is interaction as well as if eventType is escape.  
            # This works.
            single.dt = -INF
            return single.dt


        # a) Try interaction.
        interactionHorizon = single.getMinRadius() * INTERACTION_HORIZON_FACTOR
        distance, surface = \
            self.getClosestSurfaceWithinRadius( single.pos, interactionHorizon,
                                                ignore=[ single, ] )
        if surface:
            obj = self.tryInteraction( single, surface )
            if obj:
                single.dt = -INF # Remove by rescheduling to past.
                return single.dt


        # b) Try pair.
        bursted = []
        burstRadius = single.getMinRadius() * SPHERE_BURST_RADIUS_FACTOR
        # Todo. Maybe make burstRadius dimension (3D/2D/1D) dependent. 
        if isinstance( single.shellList[0], Cylinder ):
            # Todo. Implement bursting of cylindrical volume using 
            # CYLINDER_BURST_RADIUS_FACTOR and CYLINDER_BURST_SIZE_FACTOR.
            # And proper collision detection that can be used by the 
            # consistency check methods also.
            # This is quick fix to not get a minimal cylinder to overlap
            # with a nearby sphere.

            burstRadius *= math.sqrt( 2 )
        neighbors = self.getNeighborsWithinRadiusNoSort( single.pos, 
                                                         burstRadius, 
                                                         ignore=[ single, ] )
        if neighbors:
            obj, bursted = self.tryPair( single, neighbors )
            if obj:
                single.dt = -INF # Remove by rescheduling to past.
                return single.dt


        # c) Try single.
        # We already know there are no other objects within getMinSize()
        # (particle radius), because propagateSingle checked it. Now check if
        # there are any within minShell, and also get closest outside of
        # minShell.
        minShell = single.getMinRadius() * MINIMAL_SINGLE_RADIUS_FACTOR,
        neighbors, _, closestObject, closestShellDistance = \
            self.getNeighbors( single.pos, minShell, ignore=[ single, ] )
        if not neighbors:
            # All neighbors are far enough away.
            self.updateSingle( single, closestObject, closestShellDistance )

            # Todo. Still need the think about why this is needed.
            bursted = uniq( bursted )
            burstedSingles = [ s for s in bursted if isinstance( s, Single ) ]
            self.restoreSingleShells( burstedSingles )

            return single.dt


        # d) Form multi. Fall back.
        obj = self.createMulti( single, neighbors )
        assert obj, 'formMulti should always succeed.'


        # Remove single by rescheduling to past.
        single.dt = -INF

        # Note for future reference. Always return single.dt for the *current* 
        # single. So never remove the *current* single from the scheduler 
        # using removeEvent.
        # PersistentIDPolicy::getIndex(): Key not  found.).
        # Also don't do stuff like single = self.propagateSingle( single ).
        return single.dt


    ##########################################################################
    # Methods for Interactions.

    def createInteraction( self, single, surface, origin, radius, 
                           orientationVector, size, particleOffset, 
                           projectedPoint, sizeOfDomain ):
        """Create interaction.

        """
        assert single.dt == 0.0 and single.getMobilityRadius() == 0.0

        particle = single.particle

        reactionTypes = self.getReactionType1( particle.species )
        interactionType = self.getInteractionType( particle.species, surface )

        interaction = surface.defaultInteractionSingle( particle, surface, 
            reactionTypes, interactionType, origin, radius, orientationVector, 
            size, particleOffset, projectedPoint, sizeOfDomain )
        interaction.initialize( self.t )

        self.addToShellMatrix( interaction )
        interaction.dt, interaction.eventType, interaction.activeDomain = \
            interaction.determineNextEvent( )
        self.addSingleEvent( interaction )

        self.removeFromShellMatrix( single )

        cylinder = interaction.shellList[0]
        log.info( '    New %s.\n'
                  '        radius = %.3g. size = %.3g. dt = %.3g.' %
                  ( interaction, cylinder.radius, cylinder.size, 
                    interaction.dt ) )

        assert self.checkObj( interaction )

        return interaction


    def tryInteraction( self, single, surface ):
        """Todo: bursting of cylindrical volume using 
        CYLINDER_BURST_RADIUS_FACTOR and CYLINDER_BURST_SIZE_FACTOR.

        """
        particle = single.particle

        # Cyclic transpose needed when calling surface.projectedPoint!
        posTransposed = cyclicTranspose( single.pos, surface.origin, 
                                         self.worldSize )
        projectedPoint, projectionDistance = \
            surface.projectedPoint( posTransposed )

        # For interaction with a planar surface, decide orientation. 
        orientationVector = cmp( projectionDistance, 0 ) * surface.unitZ 
        particleDistance = abs( projectionDistance )

        log.debug( '        *tryInteraction: %s +\n'
                   '                         %s. particleDistance = %.10g' %
                   ( particle, surface, particleDistance ) )

        particleRadius = particle.species.radius


        # Initialize dr, dzLeft, dzRight to maximum allowed values.
        #
        # For an interaction with a PlanarSurface (membrane):
        # * dr is the radius of the cylinder.
        # * dzLeft determines how much the cylinder is sticking out on the 
        # other side of the surface, measured from the projected point.
        # * dzRight is the distance between the particle and the edge of the 
        # cylinder in the z direction.
        #
        # For interaction with a CylindricalSurface (dna):
        # * dr is the distance between the particle and the edge of the 
        # cylinder in the r direction.
        # * dzLeft is the distance from the projected point to the left edge 
        # of the cylinder.
        # * dzRight is the distance from the projected point to the right edge 
        # of the cylinder.

        if isinstance( surface, PlanarSurface ):
            dr = self.getMaxShellSize() # Free diffusion, so like free single.
            dzLeft = self.getMaxShellSize() # Doesn't matter much.

            # Todo. Note that after an unbinding,
            # r0 = MINIMAL_SEPERATION_FACTOR * (particleRadius + surface.Lz), 
            # which is very small. For now: allow for a domain that has a size 
            # of: r0 + ( MAX_DOMAIN_SIZE_FACTOR * r0 ).
            # The proper solution is to implement Greens function where either 
            # boundary is ignored, like in the 3D case.
            r0 = particleDistance - particleRadius - surface.Lz
            dzRight = particleRadius + r0 * MAX_DOMAIN_SIZE_FACTOR

            # Make sure the cylinder stays within 1 cell.
            # Todo. Assert that other surfaces are not too close.
            dzRight = min( dzRight, 
                           self.getMaxShellSize() * 2 - particleDistance )

        elif isinstance( surface, CylindricalSurface ):
            # Todo. Assert that other surfaces are not too close.
            dr = surface.radius + particleRadius * \
                ( MINIMAL_SINGLE_RADIUS_FACTOR + SURFACE_SAFETY_ZONE ) + \
                particleRadius * MINIMAL_SEPERATION_FACTOR - particleDistance 

            dzLeft = self.getMaxShellSize() # Free diffusion, like free single.
            dzRight = self.getMaxShellSize() # Free diffusion, like free single.

        # Miedema's algorithm.
        dr, dzLeft, dzRight = \
            self.determineOptimalCylinder( single, projectedPoint,
                                           orientationVector, dr, dzLeft, 
                                           dzRight, surface, particleDistance )

        # Decide minimal dr, dzleft, dzright.
        if isinstance( surface, PlanarSurface ):
            # Free diffusion in r, same as for single.
            mindr = particleRadius * MINIMAL_SINGLE_RADIUS_FACTOR
            # Leave enough for particle itself.
            mindzLeft = particleRadius * SAFETY
            # Todo. Make sure that after an escape through the z-domain of a 
            # PlanarSurfaceInteraction there is at least enough space to make 
            # a spherical single. Use MINIMAL_SINGLE_RADIUS_FACTOR and 
            # SURFACE_SAFETY_ZONE.
            mindzRight = particleRadius * MINIMAL_SEPERATION_FACTOR

        elif isinstance( surface, CylindricalSurface ):
            # Todo. Make sure that after an escape through the r-domain of a 
            # CylindricalSurfaceInteraction there is at least enough space to 
            # make a spherical single. Use MINIMAL_SINGLE_RADIUS_FACTOR and 
            # SURFACE_SAFETY_ZONE.
            mindr = particleRadius * MINIMAL_SEPERATION_FACTOR
            # Free diffusion in z, same as for single.
            mindzLeft = particleRadius * MINIMAL_SINGLE_RADIUS_FACTOR
            mindzRight = particleRadius * MINIMAL_SINGLE_RADIUS_FACTOR

        if dr < mindr or dzLeft < mindzLeft or dzRight < mindzRight:
            log.debug( '        *Interaction not possible:\n'
                       '            %s +\n'
                       '            %s.\n'
                       '            dr = %.3g. mindr = %.3g.\n'
                       '            dzLeft = %.3g. mindzLeft = %.3g.\n'
                       '            dzRight = %.3g. mindzRight = %.3g.' %
                       ( single, surface, dr, mindr, dzLeft, mindzLeft, 
                         dzRight, mindzRight ) )
            return None

        # Compute radius and size of new cylinder.

        # Todo. Let radius and size of cylinder also depend on distance 
        # towards the closest particle, taking into account the diffusion 
        # constants like in updateSingle.
        if isinstance( surface, PlanarSurface ):
            radius = dr
            # On the other side (left side) than the particle's side, make 
            # sure there is just enough space for the particle to stick out, 
            # but not more.
            dzLeft = mindzLeft 
            # Heads up: sizeOfDomain is really different from size of 
            # cylinder!
            sizeOfDomain = particleDistance + dzRight - surface.Lz
            size = ( particleDistance + dzLeft + dzRight ) / 2
        elif isinstance( surface, CylindricalSurface ):
            radius = dr + particleDistance
            # sizeOfDomain is size of cylinder * 2.
            sizeOfDomain = dzLeft + dzRight
            size = ( sizeOfDomain ) / 2

        # Compute new origin of cylinder.
        shiftZ = size - dzLeft
        origin = projectedPoint + shiftZ * orientationVector

        if isinstance( surface, PlanarSurface ):
            # Compute new particle offset relative to origin of the domain in 
            # the z-direction (z = 0), which is *at* the surface boundary. (z 
            # = L) is at the end of the cylinder.
            # (minRadius correction in the constructor).
            particleOffset = [ 0, particleDistance - surface.Lz ]
        if isinstance( surface, CylindricalSurface ):
            # Compute new particle offset in r and z-direction relative to 
            # origin of new cylinder. In the z-direction, the origin (z = 0) 
            # is at the left boundary. z = L is at the right boundary.
            # (minRadius correction in the constructor).
            particleOffset = [ particleDistance, dzLeft ]

        log.debug( '        *createInteraction\n'
                   '            dzLeft = %.3g. dzRight = %.3g.\n'
                   '            particleOffset = [ %.3g, %.3g ].\n'
                   '            sizeOfDomain = %.3g.' % 
                   ( dzLeft, dzRight, particleOffset[0], particleOffset[1], 
                     sizeOfDomain ) )

        return self.createInteraction( single, surface, origin, radius,
                                       orientationVector, size, particleOffset, 
                                       projectedPoint, sizeOfDomain )


    def determineOptimalCylinder( self, single, projectedPoint, 
                                  orientationVector, dr, dzLeft, dzRight, 
                                  surface=None, particleDistance=None ):
        """Find optimal cylinder around particle and surface, such that it 
        is not interfering with other shells. Miedema's algorithm.

        """
        allNeighbors = \
            self.getNeighborsWithinRadiusNoSort( projectedPoint, 
                                                 self.getMaxShellSize(),
                                                 ignore=[ single, ] )

        for object in allNeighbors:
            for shell in object.shellList:
                # Also works if shell is a cylinder around the same 
                # CylindricalSurface, because getRadius magically returns the 
                # size of the cylinder then.
                objectRadius = shell.radius

                # Make bursted singles look bigger, like in formPair, because 
                # the size of their shell is only particle.radius (not yet 
                # multiplied by MINIMAL_SINGLE_RADIUS_FACTOR) (and no we can 
                # not do that immediately after they are bursted, singles 
                # might start overlapping).
                if object.dt == 0.0 and object.getD() > 0:
                    # This is one of the bursted singles.
                    # Or a particle that just escaped it's multi.
                    objectRadius *= MINIMAL_SINGLE_RADIUS_FACTOR
                    #assert bursted.__contains__( object ), 'bursted = %s does 
                    #not contain %s. radius = %.3g.' % (bursted, object, 
                    #object.radius)

                objectVector = shell.origin - projectedPoint

                # Calculate zi and ri for this object.
                zi = numpy.dot( objectVector, orientationVector )
                vectorZ = zi * numpy.array( orientationVector )
                vectorR = numpy.array( objectVector ) - numpy.array( vectorZ )
                ri = numpy.linalg.norm( vectorR )

                # Calculate dr_i for this object.
                dr_i = ri - objectRadius
                if isinstance( surface, CylindricalSurface ):
                    # Only when this cylinder will be used for an interaction 
                    # with a CylindricalSurface.
                    dr_i -= particleDistance

                # Calculate dzLeft_i or dzRight_i (both are usually positive 
                # values).
                if zi < 0:
                    # Calculate dzLeft_i for this object.
                    dzLeft_i = - zi - objectRadius

                    # Miedema's algorithm left side.
                    if dzLeft_i < dzLeft and dr_i < dr:
                        if dzLeft_i > dr_i:
                            dzLeft = dzLeft_i
                        else:
                            dr = dr_i
                else:
                    # Calculate dzRight_i for this object.
                    dzRight_i = zi - objectRadius

                    if isinstance( surface, PlanarSurface ):
                        # On the particle side (right side), do Miedema's 
                        # algorithm relative to the particle's position in the 
                        # z direction.
                        # Only when this cylinder will be used for an 
                        # interaction with a PlanarSurface.
                        dzRight_i -= particleDistance

                    # Miedema's algorithm right side.
                    if dzRight_i < dzRight and dr_i < dr:
                        if dzRight_i > dr_i:
                            dzRight = dzRight_i
                        else:
                            dr = dr_i

        return dr, dzLeft, dzRight


    def burstSingle( self, single ):
        assert self.t >= single.lastTime
        assert self.t <= single.lastTime + single.dt
        assert single.radius >= single.getMinRadius()

        # Make sure propagateSingle thinks this is an escape. Only needed for 
        # interactionSingles because of the way propagateSingle works. Maybe 
        # make nicer later.
        single.eventType = EventType.ESCAPE
        newsingle = self.propagateSingle( single )

        if isinstance( single, InteractionSingle ):
            # Removing the event has to be done for *bursting* 
            # *Interaction*Singles, not for propagating InteractionSingles nor 
            # for bursting NonInteractionSingles.
            self.removeEvent( single )
        else:
            self.updateEvent( self.t, single )

        return newsingle


    def propagateSingle( self, single ):
        """The difference between a burst and a propagate is that a burst 
        always takes place before the actual scheduled event for the single, 
        while propagateSingle can be called for an escape event.

        Another subtle difference is that burstSingle always reschedules 
        (updateEvent) the single, while just calling propagate does not. This 
        works since single.dt should be returned to the scheduler after 
        calling propagateSingle.

        The return value is only used if called from burstSingle.

        """
        newpos = single.drawNewPosition( self.t - single.lastTime ) 
        self.applyBoundary( newpos )
        if not self.checkOverlap( newpos, single.getMinRadius(),
                                  ignore=[ single.particle ] ):
            raise RuntimeError( 'propagateSingle: checkOverlap failed.' )

        self.moveParticle( single.particle, newpos )

        if single.eventType == EventType.ESCAPE:
            if isinstance( single, InteractionSingle ):
                # For reactions and interactions we create a new single and 
                # get rid of the old interactionSingle in fireSingleReaction. 
                # For escapes (and bursts, because burstSingle sets 
                # single.eventType to ESCAPE) we do it here.
                self.removeFromShellMatrix( single )
                newsingle = self.createSingle( single.particle )
                log.info( '    *New %s.\n'
                          '        radius = %.3g. dt = %.3g.' %
                          ( newsingle, newsingle.radius, newsingle.dt ) )
                return newsingle
            else:
                single.pos = newpos
                single.initialize( self.t )
                self.updateShellMatrix( single )
                return single
        else:
            # REACTION. No need to update, single is removed anyway.
            pass


    def restoreSingleShells( self, singles ):
        for single in singles:
            assert single.isReset()
            closestObject, distanceToShellOfClosestObject = \
                self.getClosestObj( single.pos, ignore=[ single, ] )

            self.updateSingle( single, closestObject, 
                               distanceToShellOfClosestObject )
            self.updateEvent( self.t + single.dt, single )
            #log.debug( '        *Restore shell %s %.3g dt = %.3g.\n'
            #           '            closest %s %.3g' %
            #           ( single, single.radius, single.dt, c, d ) )


    def updateSingle( self, single, closestObject, distanceToShell ): 
        """Draw new shell + new event time.

        """
        # InteractionSingles (cylinders) should never be updated.
        assert not isinstance( single, InteractionSingle )

        if isinstance( closestObject, Single ):
            shellSize = self.calculateSingleShellSize( single, closestObject, 
                                                       distanceToShell )
        else:
            # closestObject != Single. Pair or Multi or Surface.
            shellSize = distanceToShell / SAFETY
            shellSize = max( shellSize, single.getMinRadius() )


        single.radius = min( shellSize, self.getMaxShellSize() )

        single.dt, single.eventType, single.activeDomain = \
            single.determineNextEvent( )
        single.lastTime = self.t

        # No need for self.updateEvent(), single.dt is returned from 
        # fireSingle(), or this is done in restoreSingleShells.
        self.updateShellMatrix( single )

        log.info( '    *updateSingle\n'
                  '        %s.\n'
                  '        radius = %.3g. dt = %.3g.\n'
                  '        closest = \n'
                  '        %s.\n'
                  '        distanceToShell = %.3g' %
                  ( single, single.radius, single.dt, closestObject, 
                    distanceToShell ) )


    def calculateSingleShellSize( self, single, single2, 
                                  distanceToShellOfSingle2 ):
        """Decide on a new shellSize for a single when the closestObj is also 
        a single. We already know that that single's shell is at least far 
        enough away to build a shell using MINIMAL_SINGLE_RADIUS_FACTOR (we 
        checked for that in fireSingle).

        Input: 
        * single to determine the new shellSize of.
        * closest Single
        * distance from here to *shell* of that Single. (which was also used 
          to actually determine the closest single).

        Action:
        Basically make the shellSize maximally 1/2 the distance between the 
        particles (depending on the diffusion constants).

        Possible improvement: Find really the closest particle, not the 
        particle that is in the closest shell, and compare to that also.

        """
        # Determine distance from here to *pos* of single2, so the distance 
        # between the particles.
        distanceBetweenParticles = self.distance( single.pos, single2.pos )

        minSize1 = single.getMinRadius()
        D1 = single.getD()
        if D1 == 0:
            return minSize1

        D2 = single2.getD()
        minSize2 = single2.getMinRadius()
        minSize12 = minSize1 + minSize2
        sqrtD1 = math.sqrt( D1 )
        shellSize = min( sqrtD1 / ( sqrtD1 + math.sqrt( D2 ) ) *
                         ( distanceBetweenParticles - minSize12 ) + minSize1, 
                         distanceToShellOfSingle2 )
        shellSize /= SAFETY
        shellSize = max( shellSize, minSize1 ) # not smaller than the radius
        return shellSize


    def fireSingleReaction( self, single, interactionFlag=False ):
        self.removeFromShellMatrix( single )

        reactantSpecies = single.particle.species
        currentSurface = single.surface
        oldpos = single.particle.pos.copy()  # Why copy?

        if interactionFlag:
            reactionType = single.interactionType
        else:
            reactionType = single.drawReactionType()

        log.debug( '        *fireSingleReaction: %s.' % ( reactionType ) )

        if len( reactionType.products ) == 0:
            # Also works for surface absorption reaction.
            self.removeParticle( single.particle )
            self.lastReaction = Reaction( reactionType,
                                          [ single.particle ], [] )

        elif len( reactionType.products ) == 1:
            productSpecies = reactionType.products[0]

            if ( self.isSurfaceUnbindingReaction( reactionType ) ):
                newpos = \
                    currentSurface.randomUnbindingSite( oldpos, 
                                                        productSpecies.radius )

                # Because we do +1 at end of this method.
                self.reactionEvents -= 1
                self.unbindingEvents += 1

                # After unbinding, always clear target volume. 
                self.clearVolume( newpos, productSpecies.radius )
            else:
                if( self.isSurfaceBindingReaction( reactionType ) ):
                    newpos, _ = productSpecies.surface.projectedPoint( oldpos )

                    # Because we do +1 at end of this method.
                    self.reactionEvents -= 1
                    self.bindingEvents += 1
                else: 
                    # No change.
                    newpos = oldpos

                # After binding or normal reaction, only clear volume if new 
                # radius is bigger.
                if reactantSpecies.radius < productSpecies.radius:
                    self.clearVolume( newpos, productSpecies.radius )

            if not self.checkOverlap( newpos, productSpecies.radius,
                                      ignore=[ single.particle, ] ):
                log.info( '    no space for product particle.' )
                raise NoSpace()

            self.removeParticle( single.particle )
            self.applyBoundary( newpos )
            newparticle = self.createParticle( productSpecies, newpos )
            newsingle = self.createSingle( newparticle )

            # lastReaction is used just for logging?
            self.lastReaction = Reaction( reactionType, [ single.particle ], 
                                          [ newparticle ] )
            log.info( '    New %s.\n        radius = %.3g. dt = %.3g.' %
                      ( newsingle, newsingle.radius, newsingle.dt ) )

        elif len( reactionType.products ) == 2:
            productSpecies1 = reactionType.products[0]
            productSpecies2 = reactionType.products[1]

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

                vector = \
                    currentSurface.randomVector( particleRadius12 *
                                                 MINIMAL_SEPERATION_FACTOR )

                # Place particles according to the ratio D1:D2
                # this way, species with D == 0 doesn't move.
                # FIXME: what if D1 == D2 == 0?
                while 1:
                    newpos1 = oldpos + vector * ( D1 / D12 )
                    newpos2 = oldpos - vector * ( D2 / D12 )
                    self.applyBoundary( newpos1 )
                    self.applyBoundary( newpos2 )

                    if self.distance( newpos1, newpos2 ) >= particleRadius12:
                        break

                    vector *= MINIMAL_SEPERATION_FACTOR
                # accept the new positions if there is enough space.
                if ( self.checkOverlap( newpos1, particleRadius1,
                                        ignore=[ single.particle, ] ) and
                     self.checkOverlap( newpos2, particleRadius2,
                                        ignore=[ single.particle, ] ) ):
                    break
            else:
                log.info( '    no space for product particles.' )
                raise NoSpace()

            self.removeParticle( single.particle )
            particle1 = self.createParticle( productSpecies1, newpos1 )
            particle2 = self.createParticle( productSpecies2, newpos2 )
            newsingle1 = self.createSingle( particle1 )
            newsingle2 = self.createSingle( particle2 )

            self.lastReaction = Reaction( reactionType, [ single.particle ], 
                                          [ particle1, particle2 ] )

            log.info( '    New %s.\n        radius = %.3g. dt = %.3g.' %
                      ( newsingle1, single.radius, single.dt ) )
            log.info( '    New %s.\n        radius = %.3g. dt = %.3g.' %
                      ( newsingle2, single.radius, single.dt ) )

        else:
            raise RuntimeError( 'num products >= 3 not supported.' )

        self.reactionEvents += 1


    ##########################################################################
    # Methods for Pairs.

    def createPair( self, single1, single2, shellSize ):
        assert single1.dt == 0.0 and single1.getMobilityRadius() == 0.0
        assert single2.dt == 0.0 and single2.getMobilityRadius() == 0.0
        assert single1.surface == single2.surface
        reactionType = self.reactionTypeMap2.get( ( single1.particle.species, 
                                                    single2.particle.species ) )
        pair = single1.surface.defaultPair( single1, single2, shellSize,
                                            reactionType, 
                                            Delegate( self, 
                                                      EGFRDSimulator.distance ),
                                            self.getWorldSize() )
        pair.initialize( self.t )

        self.removeFromShellMatrix( single1 )
        self.removeFromShellMatrix( single2 )
        self.addToShellMatrix( pair )

        # Formerly known as the impure function Pair.determineNextEvent().
        dtSingleReaction, reactingSingle = pair.drawSingleReactionTime( )
        dtEscape, activeDomain = pair.drawEscapeOrReactionTime( )
        assert type(dtSingleReaction).__name__ == 'float'
        assert type(dtEscape).__name__ == 'float'

        if dtSingleReaction < dtEscape:
            # This stands for *single* reaction.
            pair.eventType = EventType.REACTION
            pair.reactingSingle = reactingSingle
            pair.dt = dtSingleReaction
        else:
            # This stands for escape or normal reaction.
            pair.eventType = EventType.ESCAPE
            pair.activeDomain = activeDomain
            pair.dt = dtEscape

        pair.lastTime = self.t
        assert pair.dt >= 0

        self.addPairEvent( pair )
        # After returning, single1 is scheduled to past in fireSingle.
        self.removeEvent( single2 )

        assert self.checkObj( pair )

        return pair


    def firePair( self, pair ):
        assert self.checkObj( pair )

        particle1 = pair.single1.particle
        particle2 = pair.single2.particle

        # Two cases:
        #  0. Single reaction
        #  1. All other cases
        # First handle *single* reaction case.
        if pair.eventType == EventType.REACTION:
            self.pair_steps[3] += 1
            reactingsingle = pair.reactingSingle
            log.info( '    Pair SINGLE REACTION' )
            if reactingsingle == pair.single1:
                theothersingle = pair.single2
            else:
                theothersingle = pair.single1

            self.burstPair( pair )
            self.addSingleEvent( theothersingle )
            log.info( '    New %s.\n        radius = %.3g. dt = %.3g.' %
                      ( theothersingle, theothersingle.radius, 
                        theothersingle.dt ) )

            try:
                self.fireSingleReaction( reactingsingle )
            except NoSpace:
                self.addToShellMatrix( reactingsingle )
                self.rejectedMoves += 1
                reactingsingle.dt = 0
                self.addSingleEvent( reactingsingle )

            pair.dt = -INF
            return pair.dt

        # All other cases
        #
        # Decide if this is a *pair* reaction (0) or an escape (1) through 
        # either r or R. 
        #
        # When this pair was initialized, pair.determineNextEvent() was 
        # called, and that set the pair.activeDomain we use here. Two 
        # possibilities. If active:
        # * IV domain:
        #   drawEventType returns either:
        #    - Reaction. Escape flag is not set (wouldn't be used anyway, 
        #    because only CoM is updated)
        #    - Escape. drawEventType sets escape flag in IV domain, fixes exit 
        #    point on IV shell, used via propagatePair.
        # * CoM domain:
        #    drawEventType returns:
        #    - Escape. drawEventType sets escape flag in CoM domain, fixes 
        #    exit point on CoM shell, used via propagatePair.
        eventType = pair.activeDomain.drawEventType( pair.dt )
        log.info( '    Pair ' + str( eventType ) )


        if eventType == EventType.REACTION:
            self.pair_steps[0] += 1
            if len( pair.rt.products ) == 1:
                species3 = pair.rt.products[0]

                # calculate new R
                newCoM = pair.drawNewCoM( pair.dt )

                assert self.distance( pair.CoM, newCoM ) + species3.radius < \
                       pair.shellSize

                self.applyBoundary( newCoM )

                self.removeParticle( particle1 )
                self.removeParticle( particle2 )

                particle = self.createParticle( species3, newCoM )
                newsingle = self.createSingle( particle )

                self.reactionEvents += 1
                self.lastReaction = Reaction( pair.rt, [ particle1, particle2 ],
                                              [ particle ] )

                log.info( '    New %s.\n        radius = %.3g. dt = %.3g.' %
                          ( newsingle, newsingle.radius, newsingle.dt ) )

            else:
                raise NotImplementedError( 'num products >= 2 not supported.' )

            self.removeFromShellMatrix( pair )
            pair.dt = -INF
            return pair.dt

        # Escaping through a_r or escaping through a_R. Make use of escape flag 
        # magic.
        if eventType == EventType.ESCAPE:
            if pair.activeDomain == pair.domains[0]:
                # CoM
                self.pair_steps[2] += 1
            else:
                # IV
                self.pair_steps[1] += 1
            self.propagatePair( pair )
        else:
            raise SystemError( 'Bug: invalid eventType.' )

        # This has to be done before the following clearVolume().
        self.removeFromShellMatrix( pair )

        # The singles were there all along, now put them in the scheduler again.
        self.addSingleEvent( pair.single1 )
        self.addSingleEvent( pair.single2 )

        assert self.checkObj( pair.single1 )
        assert self.checkObj( pair.single2 )

        pair.dt = -INF
        return pair.dt


    def burstPair( self, pair ):
        assert self.t >= pair.lastTime
        assert self.t <= pair.lastTime + pair.dt

        single1, single2 = pair.single1, pair.single2

        if self.t - pair.lastTime > 0.0:
            self.propagatePair( pair )
        else:
            single1.initialize( self.t )
            single2.initialize( self.t )
            self.addToShellMatrix( single1 )
            self.addToShellMatrix( single2 )

        self.removeFromShellMatrix( pair )
        return single1, single2


    def propagatePair( self, pair ):
        single1, single2 = pair.single1, pair.single2
        single1.pos, single2.pos = pair.drawNewPositions( self.t - 
                                                          pair.lastTime )
        self.applyBoundary( single1.pos )
        self.applyBoundary( single2.pos )
        assert self.checkOverlap( single1.pos, single1.getMinRadius(),
                                  ignore=[ single1.particle, 
                                           single2.particle ] )
        assert self.checkOverlap( single2.pos, pair.single2.getMinRadius(),
                                  ignore=[ single1.particle,
                                           single2.particle ] )
        assert pair.checkNewpos( single1.pos, single2.pos )
        self.moveParticle( single1.particle, single1.pos )
        self.moveParticle( single2.particle, single2.pos )
        single1.initialize( self.t )
        single2.initialize( self.t )
        self.addToShellMatrix( single1 ) # Add, not update.
        self.addToShellMatrix( single2 )
        log.info( '    New %s.\n        radius = %.3g. dt = %.3g.' %
                  ( single1, single1.radius, single1.dt ) )
        log.info( '    New %s.\n        radius = %.3g. dt = %.3g.' %
                  ( single2, single2.radius, single2.dt ) )


    def tryPair( self, single1, neighbors ):
        """Decide if pair makes sense, and create it if so.

        """
        assert neighbors
        bursted = self.burstNonMultis( neighbors )

        # Assume the first neighbor is is also the closest neighbor, although 
        # after burstNonMultis we can not be sure.
        single2 = bursted.pop(0)

        if isinstance( single2, Multi ):
            # Don't try to make a pair with a multi.
            return None, bursted
        elif isinstance( single2, Single ):
            # Try forming a Pair only if singles are on same surface.
            if single1.surface != single2.surface:
                log.debug( '        *Pair( %s,\n'
                           '               %s )\n'
                           '            not formed: not on same surface.' %
                           ( single1.particle, single2.particle ) )
                return None, bursted
        else:
            assert False, 'Do not reach here.'


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
        distanceFromSigma = pairDistance - sigma
        assert distanceFromSigma >= 0, \
               '(pair gap) between %s and %s = %.3g < 0' % \
               ( single1, single2, distanceFromSigma )

        shellSize1 = pairDistance * D1 / D12 + radius1
        shellSize2 = pairDistance * D2 / D12 + radius2
        shellSizeMargin1 = radius1 * 2 #* SINGLE_SHELL_FACTOR
        shellSizeMargin2 = radius2 * 2 #* SINGLE_SHELL_FACTOR
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
                            distanceFromSigma * 100 + sigma + shellSizeMargin )

        if minShellSizeWithMargin >= maxShellSize:
            log.debug( '        *Pair( %s,\n'
                       '               %s )\n'
                       '            not formed: minShellSize >= maxShellSize.' %
                       ( single1.particle, single2.particle ) )
            return None, bursted

        # Here, we have to take into account of the bursted Singles in this 
        # step.  The simple check for closest below could miss some of them, 
        # because sizes of these Singles for this distance check has to 
        # include MINIMAL_SINGLE_RADIUS_FACTOR, while these bursted objects 
        # have zero mobility radii. This is not beautiful, a cleaner framework 
        # may be possible.

        closest, closestShellDistance = DummySingle(), INF
        for b in bursted:
            if isinstance( b, Single ):
                d = self.distance( com, b.pos ) - \
                    b.getMinRadius() * MINIMAL_SINGLE_RADIUS_FACTOR
                if d < closestShellDistance:
                    closest, closestShellDistance = b, d

        if closestShellDistance <= minShellSizeWithMargin:
            log.debug( '        *Pair( %s,\n'
                       '               %s )\n'
                       '            not formed: squeezed by bursted neighbor\n'
                       '            %s.' %
                       ( single1.particle, single2.particle, closest ) )
            return None, bursted

        c, d = self.getClosestObj( com, ignore=[ single1, single2 ] )
        if d < closestShellDistance:
            closest, closestShellDistance = c, d

        log.debug( '        *Try pair. Closest neighbor =\n'
                   '            %s.\n'
                   '            closestShellDistance = %.3g.\n'
                   '            minShellSize = %.3g.' %
                   ( closest, closestShellDistance, minShellSizeWithMargin ) )

        if isinstance( closest, Single ):

            D_closest = closest.particle.species.D
            D_tot = D_closest + D12
            closestDistance = self.distance( com, closest.pos )

            closestMinSize = closest.getMinRadius()
            closestMinShell = closestMinSize * MINIMAL_SINGLE_RADIUS_FACTOR

            shellSize = min( ( D12 / D_tot ) *
                             ( closestDistance - minShellSize 
                               - closestMinSize ) + minShellSize,
                             closestDistance - closestMinShell,
                             closestShellDistance )

            shellSize /= SAFETY
            assert shellSize < closestShellDistance

        else:
            assert isinstance( closest, ( Pair, Multi, Surface, DummySingle ) ), \
                   str( closest )

            shellSize = closestShellDistance / SAFETY

        if shellSize <= minShellSizeWithMargin:
            log.debug( '        *Pair( %s,\n'
                       '               %s )\n'
                       '            not formed: squeezed by\n'
                       '            %s.\n'
                       '            shellSize = %.3g. minShellSize = %.3g.' %
                       ( single1.particle, single2.particle, closest, 
                         shellSize, minShellSizeWithMargin ) )
            return None, bursted

        d1 = self.distance( com, single1.pos )
        d2 = self.distance( com, single2.pos )

        if shellSize < 1.3 * \
                       max( d1 + single1.getMinRadius() * 
                                 MINIMAL_SINGLE_RADIUS_FACTOR,
                            d2 + single2.getMinRadius() * 
                                 MINIMAL_SINGLE_RADIUS_FACTOR ): 
            log.debug( '        *Pair( %s,\n'
                       '               %s )\n'
                       '            not formed: singles are better.' % 
                       ( single1.particle, single2.particle ) )
            return None, bursted

        shellSize = min( shellSize, maxShellSize )

        assert closestShellDistance == INF or shellSize < \
               closestShellDistance
        assert shellSize >= minShellSizeWithMargin
        assert shellSize <= maxShellSize

        # 3. Ok, Pair makes sense.  Create one.
        pair = self.createPair( single1, single2, shellSize )

        log.info( '    New %s.\n        radius = %.3g. pairDistance = %.3g. '
                  'dt = %.3g.' % 
                  ( pair, pair.shellSize, pairDistance, pair.dt ) )

        return pair, bursted


    ##########################################################################
    # Methods for Multis.

    def fireMulti( self, multi ):
        sim = multi.sim
        self.multi_steps[2] += 1  # multi_steps[2]: total multi steps
        sim.step()
        #sim.sync()

        if sim.lastReaction:
            log.info( '    Bd reaction.' )

            self.burstMulti( multi )
            self.reactionEvents += 1
            self.lastReaction = sim.lastReaction
            self.multi_steps[EventType.REACTION] += 1
            return -INF

        if sim.escaped:
            log.info( '    Multi particle escaped.' )

            self.burstMulti( multi )
            self.multi_steps[EventType.ESCAPE] += 1
            return -INF

        #log.info( '    multi stepped %d steps, duration %.3g, dt = %.3g' %
        #          ( additionalSteps + 1, sim.t - startT + sim.dt, dt ) )
        return multi.dt


    def burstMulti( self, multi ):
        self.removeFromShellMatrix( multi )
        singles = []
        for particle in multi.sim.particleList:
            single = self.createSingle( particle )
            singles.append( single )
        return singles


    def createMulti( self, single, neighbors ):
        neighborsString = ''
        for n in neighbors:
            neighborsString += str( n ) + ',\n            '
        log.debug( '        *Try to form Multi: %s +\n'
                   '            [ %s ].' % ( single, neighborsString ) )

        closest = neighbors[0]
        if isinstance( closest, Single ) or isinstance( closest, Surface ):

            multi = Multi( self )
            #multi.initialize( self.t )

            # Use addToMulti here instead of addToMultiRecursive because we 
            # already found this single's neighbors (we don't want to do that 
            # twice).
            self.addToMulti( single, multi )
            self.removeFromShellMatrix( single )
            for neighbor in neighbors:
                self.addToMultiRecursive( neighbor, multi )

            multi.initialize( self.t )

            self.addToShellMatrix( multi )
            self.addMultiEvent( multi )
            return multi#, bursted

        elif isinstance( closest, Multi ):

            multi = closest
            log.info( '    multi merge %s,\n        %s' % ( single, multi ) )

            self.removeFromShellMatrix( multi )

            self.addToMulti( single, multi )
            self.removeFromShellMatrix( single )
            for neighbor in neighbors[1:]:
                self.addToMultiRecursive( neighbor, multi )

            multi.initialize( self.t )

            # Add 1 shell for each particle to main simulator's shellMatrix.
            self.addToShellMatrix( multi )
            self.updateEvent( self.t + multi.dt, multi )
            return multi#, bursted

        assert False, 'do not reach here'


    def addToMultiRecursive( self, obj, multi ):
        """Add obj to multi using addToMulti().

        If obj is a Single, also add any neighbors that lie within 
        shellSize of obj, where shellSize is MULTI_SHELL_FACTOR times the 
        radius of the particle.

        If obj is itself a Multi, merge them.

        """
        if isinstance( obj, Surface ):
            # No need to add other particles recursively, surface doesn't 
            # move.
            multi.addSurface( obj )
        elif isinstance( obj, Single ):
            if obj.particle in multi.sim.particleList:  # Already in the Multi.
                return
            assert obj.isReset()

            # Add particle + new shell to multi.
            self.addToMulti( obj, multi )
            self.removeFromShellMatrix( obj )
            self.removeEvent( obj )

            # Find any neighbouring particles that lie within shellSize and 
            # also add them (recursively) to the multi.
            shellSize = obj.particle.species.radius * MULTI_SHELL_FACTOR
            neighbors = self.getNeighborsWithinRadiusNoSort( obj.pos, shellSize,
                                                             ignore=[ obj, ] )
            bursted = self.burstNonMultis( neighbors )
            neighborDists = self.objDistanceArray( obj.pos, bursted )
            neighbors = [ bursted[i] for i in 
                          ( neighborDists <= shellSize ).nonzero()[0] ]

            for neighbor in neighbors:
                self.addToMultiRecursive( neighbor, multi )

        elif isinstance( obj, Multi ):
            if not obj.sim.particleList[0] in multi.sim.particleList:
                self.mergeMultis( obj, multi )
                self.removeFromShellMatrix( obj )
                self.removeEvent( obj )
            else:
                log.debug( '        *%s already added. skipping.' % obj )
        else:
            assert False, 'do not reach here.'  # Pairs are bursted



    def addToMulti( self, single, multi ):
        """Add to multi object:
            - the particle
            - a shell a bit bigger than radius (using MULTI_SHELL_FACTOR)

        """
        log.info( '    Adding %s to\n        %s' % ( single, multi ) )

        shellSize = single.particle.species.radius * MULTI_SHELL_FACTOR
        multi.addParticle( single.particle )
        multi.addShell( single.pos, shellSize )


    def mergeMultis( self, multi1, multi2 ):
        """Merge multi1 into multi2.

        """
        log.info( '    merging %s to %s' % ( multi1, multi2 ) )
        assert not multi1.sim.particleList[0] in multi2.sim.particleList

        for i, particle in enumerate( multi1.sim.particleList ):
            # FIXME: shells should be renewed
            multi2.addParticle( particle )
            shell = multi1.shellList[i]
            multi2.addShell( shell.origin, shell.radius )

        multi2.initialize( self.t )



    #
    # statistics reporter
    #

    def print_report(self, out=None):
        report = '''
t = %g
steps = %d 
\tSingle:\t%d\t(escape: %d, reaction: %d)
\tPair:\t%d\t(escape r: %d, R: %d, reaction pair: %d, single: %d)
\tMulti:\t%d\t(escape: %d, reaction: %d)
total reactions = %d
rejected moves = %d
'''\
            % (self.t, self.stepCounter,
               numpy.array(self.single_steps.values()).sum(),
               self.single_steps[EventType.ESCAPE],
               self.single_steps[EventType.REACTION],
               numpy.array(self.pair_steps.values()).sum(),
               self.pair_steps[1],
               self.pair_steps[2],
               self.pair_steps[0],
               self.pair_steps[3],
               self.multi_steps[2], # total multi steps
               self.multi_steps[EventType.ESCAPE],
               self.multi_steps[EventType.REACTION],
               self.reactionEvents,
               self.rejectedMoves
               )

        print >> out, report
