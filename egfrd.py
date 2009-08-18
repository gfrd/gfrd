#!/usr/env python
import math
import numpy

from gfrdbase import *

from single import *
from pair import *
from multi import *
from shape import *

from itertools import izip
from vtklogger import VTKLogger

from log import *


'''
Returning a dt to the scheduler reschedules the event, -INF removes it.
'''
class EGFRDSimulator( ParticleSimulatorBase ):
    def __init__( self, worldSize, logdir = None ):
        if logdir:
            log.info( '\tStarted vtk logger: ' + logdir )
            self.vtklogger = VTKLogger( self, logdir )
        else:
            self.vtklogger = None

        self.sphereMatrix = SphereMatrix()
        self.cylinderMatrix = CylinderMatrix()
        self.boxMatrix = BoxMatrix() # Not used.
        self.objectMatrices = [ self.sphereMatrix, self.cylinderMatrix, self.boxMatrix ]
        #self.sm2 = pObjectMatrix()

        ParticleSimulatorBase.__init__( self, worldSize )

        self.MULTI_SHELL_FACTOR = 0.05 # Should be smaller than SINGLE_SHELL_FACTOR!
        self.SINGLE_SHELL_FACTOR = 0.1

        self.isDirty = True
        self.scheduler = EventScheduler()

        self.smallT = 1e-8  # FIXME: is this ok?
        self.userMaxShellSize = INF
        self.reset()


    # Called from ParticleSimulatorBase.__init__()
    def setWorldSize( self, size ):
        if isinstance( size, list ) or isinstance( size, tuple ):
            size = numpy.array( size )

        ParticleSimulatorBase.setWorldSize( self, size )

        self.sphereMatrix.setWorldSize( size )
        self.cylinderMatrix.setWorldSize( size )
        self.boxMatrix.setWorldSize( size )


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
	self.interactionEvents = 0
        self.lastEvent = None
        self.lastReaction = None
        self.isDirty = True
        self.errors = 0
        #self.initialize()


    def initialize( self, a=None ):
        ParticleSimulatorBase.initialize( self )
        self.setAllRepulsive()
        self.scheduler.clear()
        self.sphereMatrix.clear()
        self.cylinderMatrix.clear()
        self.boxMatrix.clear()

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
            self.initialize(1)
        if self.vtklogger:
            self.vtklogger.log()

        self.lastReaction = None
            
        #if self.stepCounter % 100 == 0:
        #    self.check()
        
        self.stepCounter += 1
        event = self.scheduler.getTopEvent()
        self.t, self.lastEvent = event.getTime(), event.getArg()

        log.info( '\n\n\t%d: t=%g dt=%g\n\tevent=%s. reactions=%d interactions=%d rejectedmoves=%d errors=%d.' 
                      % ( self.stepCounter, self.t, self.dt, self.lastEvent, 
                          self.reactionEvents, self.interactionEvents, self.rejectedMoves, self.errors ) )
        
        self.scheduler.step()

        nextTime = self.scheduler.getTopTime()
        self.dt = nextTime - self.t

        # assert if not too many successive dt=0 steps occur.
        if self.dt == 0.0:
            self.zeroSteps += 1
            if self.zeroSteps >= max( self.scheduler.getSize() * 3, 10 ):
                raise Stop( 'Too many dt=zero steps.  simulator halted?' )
        else:
            self.zeroSteps = 0
        assert self.scheduler.getSize() != 0


    def stop( self, t ):
        if self.vtklogger:
            self.vtklogger.stop()
        log.info( '\tStop at %g.' % t )

        if feq( self.t, t ):
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
                log.debug( '\t\tDebug. burst %s, lastTime= %g' % 
                           ( str( obj ), obj.lastTime ) )
                self.burstSingle( obj )
            else:
                assert False, 'do not reach here'


        # then burst all Pairs and Multis.
        log.debug( '\t\tDebug. burst %s' % nonSingleList )
        self.burstObjs( nonSingleList )
        self.dt = 0.0

    
    ##########################################################################
    '''
    Burst methods.
    '''
    def clearVolume( self, pos, radius, ignore=[] ):
        neighbors = self.getNeighborsWithinRadiusNoSort( pos, radius, ignore )
        return self.burstObjs( neighbors )


    def burstNonMultis( self, neighbors ):
        bursted = []
        for obj in neighbors:
            if not (isinstance( obj, Multi ) or isinstance( obj, Surface )):
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
        log.info( '\tBursting %s ' % ( str( obj ) ) )
        if isinstance( obj, Single ):
            obj = self.burstSingle( obj )
            return [obj,]
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
            raise NotImplementedError


    ##########################################################################
    '''
    Event scheduler stuff.
    '''
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
    '''
    Shell matrix stuff. 
    '''
    def setMatrixSize( self, size ):
        ParticleSimulatorBase.setMatrixSize( self, size )
        self.sphereMatrix.setMatrixSize( size )
        self.cylinderMatrix.setMatrixSize( size )
        self.boxMatrix.setMatrixSize( size )


    def getMatrixCellSize( self ):
        assert self.sphereMatrix.cellSize == self.cylinderMatrix.cellSize == self.boxMatrix.cellSize
        return self.sphereMatrix.cellSize


    '''
    Add all shells from obj.shellList to one of the shell matrices.
    '''
    def addToShellMatrix( self, obj ):
        for i, shell in enumerate( obj.shellList ):
            key = (obj, i)
            if isinstance( shell, Sphere ):
                self.sphereMatrix.add( key, shell.origin, shell.radius )
            elif isinstance( shell, Cylinder ):
                self.cylinderMatrix.add( key, shell )
            elif isinstance( shell, Box ):
                self.boxMatrix.add( key, shell )
            else: raise KeyError, 'Objecttype does not exit'


    def updateShellMatrix( self, obj ):
        for i, shell in enumerate( obj.shellList ):
            key = (obj, i)
            if isinstance( shell, Sphere ):
                self.sphereMatrix.update( key, shell.origin, shell.radius )
            elif isinstance( shell, Cylinder ):
                self.cylinderMatrix.update( key, shell )
            elif isinstance( shell, Box ):
                self.boxMatrix.update( key, shell )
            else: raise KeyError, 'Objecttype does not exit'


    def removeFromShellMatrix( self, obj ):
        for i, shell in enumerate( obj.shellList ):
            key = (obj, i)
            if isinstance( shell, Sphere ):
                self.sphereMatrix.remove( key )
            elif isinstance( shell, Cylinder ):
                self.cylinderMatrix.remove( key )
            elif isinstance( shell, Box ):
                self.boxMatrix.remove( key )
            else: raise KeyError, 'Objecttype does not exit'

    '''
    Find closest n shells.
    These methods return a tuple ( neighbors, distances ).
    Options:
        * Sort yes/no.
        * Within radius yes/no.
        * With ignore list yes/no.
    '''
    def getClosestObj( self, pos, ignore=[] ):
        # No radius (i.e. all), no ignore list, keys not extracted.
        '''
        Don't be smart and do:
            keys, distances, _, _ = self.getNeighbors( pos, ignore )
            return keys[0], distance[0]
        since that would unneccesarily iterate one more time over all 
        objects in self.getNeighbors. 
        '''
        closestObject = DummySingle()
        closestDistance = INF

        for objectMatrix in self.objectMatrices:
            keys, distances = objectMatrix.getNeighbors( pos )
            '''
            Don't be clever and think you can do:
                keys, distances = objectMatrix.getNeighbors( pos, 1 )
            Because then you are ignoring the ignore list.
            '''
            for i, key in enumerate( keys ):
                single = key[0]
                if single not in ignore and distances[i] < closestDistance:
                    closestObject, closestDistance = single, distances[i]
                    '''
                    Found yet a closer single. Break out of inner for loop 
                    and check other objectMatrices.
                    '''
                    break   

        # Surface detection.
        distanceToSurface, closestSurface = self.getClosestSurface( pos, ignore ) 
        if distanceToSurface < closestDistance:
            closestDistance = distanceToSurface
            closestObject = closestSurface

        return closestObject, closestDistance


    def getNeighborsWithinRadiusNoSort( self, pos, radius, ignore=[] ):
        neighbors = []
        for objectMatrix in self.objectMatrices:
            keys, _ =\
                objectMatrix.getNeighborsWithinRadiusNoSort( pos, radius )
            neighbors.extend( uniq( [ key[0] for key in keys if key[0] not in ignore ] ) )
        return neighbors


    '''
    Returns all singles and the distances towards their *shell* if that shell is 
    within 'radius' of 'pos', plus the closest single and the distance towards it 
    shell that is just outside of radius.
    '''
    def getNeighbors( self, pos, radius=INF, ignore=[] ):
        # Diccionaries are more efficient for lookups.
        seen = dict.fromkeys( ignore )
        neighbors = []
        distances = []

        closestSingle = DummySingle()
        closestDistance = INF

        for objectMatrix in self.objectMatrices:
            keys, dists = objectMatrix.getNeighbors( pos )

            for i, key in enumerate( keys ):
                single = key[0]
                if not single in seen:
                    '''
                    Since a single can in theory have more than 1 shell, and 
                    for each shell there is an entry in the objectMatrix, we 
                    signal here that we have already seen this single.
                    This used to be a bug: seen[ key ] = None
                    '''
                    seen[ single ] = None
                    if dists[i] > radius:
                        '''
                        This is a single that has a shell that is more than 
                        radius away from pos. If it is closer than the 
                        closest such one we found so far: store it.
                        Always break out of the inner for loop now and check 
                        the other objectMatrices.
                        '''
                        if  dists[i] < closestDistance:
                            closestSingle = single
                            closestDistance = dists[i]
                            break
                        else:
                            break # Just to be sure you get the point.
                    else:
                        neighbors.append( single )
                        distances.append( dists[i] )

        # Surface detection.
        distanceToSurface, closestSurface = self.getClosestSurface( pos, ignore ) 
        if distanceToSurface < 0:
            raise Stop( 'distanceToSurface < 0' )
        if distanceToSurface < radius:
            neighbors.append( closestSurface )
            distances.append( distanceToSurface )
        elif distanceToSurface < closestDistance:
            closestDistance = distanceToSurface
            closestSingle = closestSurface

        return neighbors, distances, closestSingle, closestDistance


    # Possible improvement: get this info from objectMatrix using 
    # objectMatrix.get().
    def objDistance( self, pos, obj ):
        if isinstance( obj, Surface ):
            posTransposed = cyclicTranspose( pos, obj.origin, self.worldSize )
            return obj.signedDistanceTo( posTransposed )
        else:
            dists = numpy.zeros( len( obj.shellList ) )
            for i, shell in enumerate( obj.shellList ):
                posTransposed = cyclicTranspose( pos, shell.origin, self.worldSize )
                dists[i] = shell.signedDistanceTo( posTransposed )
            return min( dists )


    def objDistanceArray( self, pos, objs ):
        dists = numpy.array( [ self.objDistance( pos, obj ) for obj in objs ] )
        return dists
            

    ##########################################################################
    '''
    consistency checkers
    '''
    def checkObj( self, obj ):
	#log.debug( '\t\tDebug. checkObj: %s' % (obj) )
        obj.check()

        allshells = [ ( obj, i ) for i in range( len( obj.shellList ) ) ]
        for i, shell in enumerate( obj.shellList ):
            closest, distance = self.getClosestObj( shell.origin,
                                                    ignore = [obj,] )
            radius = shell.radius


            assert radius <= self.getUserMaxShellSize(),\
                '%s shell radius larger than user-set max shell radius' % \
                str( ( obj, i ) )
            assert radius <= self.getMaxShellSize(),\
                '%s shell radius larger than simulator cell radius / 2. %g vs %g' % \
                (str( ( obj, i ) ), radius, self.getMaxShellSize())

            # Todo.
            break
            
            if not ( isinstance( obj, InteractionSingle ) and isinstance( closest, Surface ) ):
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
        if self.worldSize != self.sphereMatrix.worldSize:
            raise RuntimeError,\
                'self.worldSize != self.sphereMatrix.worldSize'

        shellPopulation = 0
        for i in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(i).getArg()
            shellPopulation += len( obj.shellList )
        
        matrixPopulation = sum(matrix.size for matrix in self.objectMatrices)
        if shellPopulation != matrixPopulation:
            raise RuntimeError,\
                'num shells %g != matrixPopulation %g' %(shellPopulation, matrixPopulation)
        
        self.sphereMatrix.check()

        for k in range( self.scheduler.getSize() ):
            obj = self.scheduler.getEventByIndex(k).getArg()
            for i in range( len( obj.shellList ) ):
                key = ( obj, i )
                if isinstance( key[0].shellList[0], Cylinder ):
                    # Todo.
                    break

                pos, radius = self.sphereMatrix.get( key )

                if ( obj.shellList[i].origin - pos ).any():
                    raise RuntimeError, \
                        '%s sphereMatrix positions consistency broken' % str( key )

                if obj.shellList[i].radius != radius:
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


    ##########################################################################
    '''
    methods for debugging.
    '''
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


    ##########################################################################
    '''
    Methods for Singles.
    '''
    def createSingle( self, particle ):
	# Poissonian reactions.
	rt = self.getReactionType1( particle.species )
	'''
	The type of defaultSingle depends on the surface this particle is on.  
	Either SphericalSingle3D, CylindricalSingle2D, or CylindricalSingle1D.
	'''
	single = particle.surface.defaultSingle( particle, rt )
	single.initialize( self.t )
	self.addToShellMatrix( single )
	self.addSingleEvent( single )
	return single


    '''
    A single particle either reacts or leaves it's shell. Here these
    events are executed.
    '''
    def fireSingle( self, single ):
	assert abs( single.dt + single.lastTime - self.t ) <= 1e-18 * self.t

        # 0. Reaction.
	if single.eventType == EventType.REACTION:
	    self.propagateSingle( single )
	    log.info( str( single.eventType ) )

	    try:
		self.fireSingleReaction( single )
	    except NoSpace:
		log.info( '\tsingle reaction; placing product failed.' )
		self.addToShellMatrix( single )
		self.rejectedMoves += 1
		single.reset()
		return single.dt

	    single.dt = -INF  # remove old single event from the Scheduler
	    return single.dt

        # Handle immobile case first.
	if single.getD() == 0:
	    # no propagation, just calculate next reaction time.
	    # (!) Not a pure function.
	    single.dt, single.eventType, single.activeDomain = single.determineNextEvent( ) 
	    single.lastTime = self.t
	    return single.dt
	
        # Propagate if dt != 0.
        if single.dt != 0.0:
            # All Singles, including NonInteractionSingles, should call 
            # drawEventType, even if it is already known that it is going to 
            # be an escape. This is because drawEventType also sets the escape 
            # flag in the active domain, which is used when calling 
            # propagateSingle and signals that in this domain the exitpoint is 
            # already known.
            single.eventType = single.activeDomain.drawEventType( single.dt )

            # Propagate this particle to the exit point on the shell.
            self.propagateSingle( single )

        if single.eventType == EventType.REACTION:
            log.info( '\tINTERACTION' )
            try:
                self.fireSingleReaction( single, True )
            except NoSpace:
                log.info( '\tsingle interaction; placing product failed.' )
                self.addToShellMatrix( single )
                self.rejectedMoves += 1
                single.reset()
                return single.dt
        else:
            log.info( '\tESCAPE' )


        if isinstance( single, InteractionSingle ):
            # A new single is created, either in propagateSingle or in 
            # fireSingleReaction. Remove the old single from scheduler.
            # For eventType is interaction as well as escapes. No other way.
            single.dt = -INF
            return single.dt


	# (2) Clear volume.
        # The single was just propagated and initialized, so it's shell has
        # size getMinSize().
	minShell = single.getMinRadius() * ( 1.0 + self.SINGLE_SHELL_FACTOR )
        # We already know there are no other objects within getMinSize()
        # (particle radius), because propagateSingle checked it. Now check if
        # there are any within minShell, and also get closest outside of
        # minShell.
	closeNeighbors, distances, closest, closestShellDistance = \
		self.getNeighbors( single.pos, minShell, ignore=[single,] )

	bursted = []
	if closeNeighbors:
            # If a closeNeighbor is already a multi, don't burst it, but
            # let it absorb the single in the next step.
	    bursted = self.burstNonMultis( closeNeighbors )
            for b in bursted:
                assert not isinstance( b, InteractionSingle )
	    obj, b = self.formInteractionOrPairOrMulti( single, bursted )
            # Now obj can be a pair or a multi, b contains all singles that 
            # were not added to it.
            # Todo. Why extend?
	    bursted.extend( b )
	    if obj:
                # Maybe restore bursted singles here? They are added to the 
                # scheduler with dt=0 already, so works this way also.
		single.dt = -INF # remove by rescheduling to past.
		return single.dt

	    # if nothing was formed, recheck closest and restore shells.
            # Maybe some other particle has come closer during the
            # bursting (while no multi had to be formed).
	    closest, closestShellDistance = \
		self.getClosestObj( single.pos, ignore = [ single, ] )

	# All neighbors are more than minShell away.
        #
        # We can get here if:
        # 1. there were closeNeighbors. But it wasn't possible to build a pair 
        # or necessary to build a multi.
        # 2. there was not a neighbor nearby.
        #
        # Update single.
	self.updateSingle( single, closest, closestShellDistance )

	bursted = uniq( bursted )
	burstedSingles = [ s for s in bursted if isinstance( s, Single ) ]
        # Probably everything would work just fine still if we didn't restore
        # here, since those burstedSingles are already in the scheduler with a
        # dt=0 and they have been given a shell with size minSize.
        #
        # Todo. Seems like this is really needed, but why?
	self.restoreSingleShells( burstedSingles )

        # Note for future reference. Always return single.dt here for the  
        # current single. So never remove the current single from the 
        # scheduler using removeEvent.
        # PersistentIDPolicy::getIndex(): Key not  found.).
        # Also don't do stuff like single = self.propagateSingle( single )
	return single.dt


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
            # We can not do this in propagate, but don't forget it.
            self.removeEvent( single )
        else:
            self.updateEvent( self.t, single )

        return newsingle

    ''' 
    The difference between a burst and a propagate is that a burst always takes 
    place before the actual scheduled event for the single, while propagate 
    can be called for an escape event.
    Another subtle difference is that burstSingle always reschedules (updateEvent) 
    the single, while just calling propagate not. This works if 
    single.dt is returned to the scheduler after calling propagateSingle.

    The return value is only used if called from burstSingle.
    '''
    def propagateSingle( self, single ):
	newpos = single.drawNewPosition(self.t - single.lastTime) 
	self.applyBoundary( newpos )
	if not self.checkOverlap( newpos, single.getMinRadius(),
                ignore = [ single.particle ] ):
            raise Stop( 'propagateSingle: checkOverlap failed.' )

	self.moveParticle( single.particle, newpos )

        if single.eventType == EventType.ESCAPE:
            if isinstance( single, InteractionSingle ):
                '''
                For reactions and interactions we create a new single and get 
                rid of the old interactionSingle in fireSingleReaction. For 
                escapes (and bursts) we do it here.
                '''
                self.removeFromShellMatrix( single )
                newsingle = self.createSingle( single.particle )
                log.info( '\tNew %s. radius=%g. dt=%g.' % ( newsingle, newsingle.radius, newsingle.dt ) )
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
	    c, d = self.getClosestObj( single.pos, ignore = [single,] )

	    self.updateSingle( single, c, d )
	    self.updateEvent( self.t + single.dt, single )
	    #log.debug( '\t\tDebug. restore shell %s %g dt %g closest %s %g' %
            #	       ( single, single.radius, single.dt, c, d ) )


    # Draw new shell + new event time.
    def updateSingle( self, single, closest, distanceToShell ): 
        # InteractionSingles (cylinders) should never have to be updated.  
        assert not isinstance( single, InteractionSingle )
	if isinstance( closest, Single ):
            # Todo.
            # and isinstance( closest.shellList[0], Sphere)
            #
            # Closest == Single.
	    distanceToClosest = self.distance( single.pos, closest.pos )
	    shellSize = self.calculateSingleShellSize( single, closest, 
						       distanceToClosest,
						       distanceToShell )
	else:
            # Closest != Single. Pair or Multi or Surface.
	    shellSize = distanceToShell / SAFETY
	    shellSize = max( shellSize, single.getMinRadius() )

	shellSize = min( shellSize, self.getMaxShellSize() )
	single.radius = shellSize

	single.dt, single.eventType, single.activeDomain = single.determineNextEvent( )
	single.lastTime = self.t
        # No need for self.updateEvent(), single.dt is returned from 
        # fireSingle(), or this is done in restoreSingleShells.
	self.updateShellMatrix( single )
	log.info( '\tUpdated %s. radius=%g. dt=%g.\n\t\tclosest=%s. distanceToShell=%s' % ( single, single.radius, single.dt, closest, distanceToShell ) )



    '''
    Decide on a new shellSize for a single when the closestObj is also a 
    single. We already know that that single's shell is at least far enough 
    away to build a shell using SINGLE_SHELL_FACTOR (we checked for that in 
    fireSingle).
    
    Input: 
    * single to determine the new shellSize of.
    * closest Single
    * distance from here to *pos* of that Single, so the distance between the 
      particles.
    * distance from here to *shell* of that Single. (which was also used to 
      actually determine the closest single).

    Action:
    Basically make the shellSize maximally 1/2 the distance between the 
    particles (depending on the diffusion constants).

    Possible improvement: Find really the closest particle, not the particle 
    that is in the closest shell, and compare to that also.
    '''
    def calculateSingleShellSize( self, single, closestSingle, 
				  distanceBetweenParticles, closestShellDistance ):
	minSize1 = single.getMinRadius()
	D1 = single.getD()
	if D1 == 0:
	    return minSize1

	D2 = closestSingle.getD()
	minSize2 = closestSingle.getMinRadius()
	minSize12 = minSize1 + minSize2
	sqrtD1 = math.sqrt( D1 )
	shellSize = min( sqrtD1 / ( sqrtD1 + math.sqrt( D2 ) )
                         * ( distanceBetweenParticles - minSize12 ) + minSize1,
			 closestShellDistance )
	shellSize /= SAFETY
	shellSize = max( shellSize, minSize1 ) # not smaller than the radius
	return shellSize


    def fireSingleReaction( self, single, interactionFlag = False ):
        self.removeFromShellMatrix( single )

	reactantSpecies = single.particle.species
	currentSurface = single.surface
	oldpos = single.particle.pos.copy()  # Why copy?
	
        if interactionFlag:
            rt = single.interactionType
        else:
            rt = single.drawReactionType()

        log.debug( '\t\tDebug. fireSingleReaction: %s.' %( rt ) )
	
	if len( rt.products ) == 0:
	    self.removeParticle( single.particle )
	    self.lastReaction = Reaction( rt, [single.particle], [] )

	elif len( rt.products ) == 1:

	    productSpecies = rt.products[0]
	    if reactantSpecies.radius < productSpecies.radius:
		self.clearVolume( oldpos, productSpecies.radius )

	    if not self.checkOverlap( oldpos, productSpecies.radius,
				      ignore = [ single.particle, ] ):
		log.info( '\tno space for product particle.' )
		raise NoSpace()

	    if isinstance( rt, SurfaceUnbindingReactionType ):
		newpos = currentSurface.randomUnbindingSite( oldpos, productSpecies.radius )
            elif isinstance( rt, SurfaceBindingInteractionType ):
                # Todo. Does this obey detailed balance?
                newpos, _ = rt.products[0].surface.projectedPoint( oldpos )

                self.reactionEvents -= 1 # Because we do +1 at end of this method.
                self.interactionEvents += 1
	    else: 
		# No change.
		newpos = oldpos

	    self.removeParticle( single.particle )
            self.applyBoundary( newpos )
	    newparticle = self.createParticle( productSpecies, newpos )
	    newsingle = self.createSingle( newparticle )

            # Todo. Is lastReaction used anywhere, or is it just logging?
	    self.lastReaction = Reaction( rt, [single.particle], [newparticle] )
            log.info( '\tNew %s. radius=%g. dt=%g.' % ( newsingle, newsingle.radius, newsingle.dt ) )

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

		vector = currentSurface.randomVector( particleRadius12 * ( 1.0 + 1e-7 ) )
	    
		'''
		place particles according to the ratio D1:D2
		this way, species with D=0 doesn't move.
		FIXME: what if D1 == D2 == 0?
		'''
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
		log.info( '\tno space for product particles.' )
		raise NoSpace()

	    self.removeParticle( single.particle )
	    particle1 = self.createParticle( productSpecies1, newpos1 )
	    particle2 = self.createParticle( productSpecies2, newpos2 )
	    newsingle1 = self.createSingle( particle1 )
	    newsingle2 = self.createSingle( particle2 )

	    self.lastReaction = Reaction( rt, [single.particle], 
					  [particle1, particle2] )

            log.info( '\tNew %s. radius=%g. dt=%g.' % ( newsingle1, single.radius, single.dt ) )
            log.info( '\tNew %s. radius=%g. dt=%g.' % ( newsingle2, single.radius, single.dt ) )

	else:
	    raise RuntimeError, 'num products >= 3 not supported.'

	self.reactionEvents += 1


    ##########################################################################
    '''
    Methods for Pairs.
    '''
    def createPair( self, single1, single2, shellSize ):
	assert single1.dt == 0.0 and single1.getMobilityRadius() == 0.0
	assert single2.dt == 0.0 and single2.getMobilityRadius() == 0.0
	assert single1.surface == single2.surface
	rt = self.reactionTypeMap2.get( (single1.particle.species, 
	    single2.particle.species) )
	pair = single1.surface.defaultPair(single1, single2, shellSize, 
		rt, Delegate( self, EGFRDSimulator.distance ), self.getWorldSize())
	pair.initialize( self.t )
	return pair


    def firePair( self, pair ):
	assert self.checkObj( pair )

	particle1 = pair.single1.particle
	particle2 = pair.single2.particle
	
	'''
	Two cases:
	  0. Single reaction
	  1. All other cases
	'''
	# First handle *single* reaction case.
	if pair.eventType == EventType.REACTION:
	    reactingsingle = pair.reactingSingle
	    log.info( '\tPair SINGLE REACTION' )
	    if reactingsingle == pair.single1:
		theothersingle = pair.single2
	    else:
		theothersingle = pair.single1

	    self.burstPair( pair )
	    self.addSingleEvent( theothersingle )
            log.info( '\tNew %s. radius=%g. dt=%g.' % ( theothersingle, theothersingle.radius, theothersingle.dt ) )

	    try:
		self.fireSingleReaction( reactingsingle )
	    except NoSpace:
		self.addToShellMatrix( reactingsingle )
		self.rejectedMoves += 1
		reactingsingle.dt = 0
		self.addSingleEvent( reactingsingle )

	    pair.dt = -INF
	    return pair.dt

	'''
	1. All other cases
	
	Decide if this is a *pair* reaction (0) or an escape (1) through either r 
	or R. 

        When this pair was initialized, pair.determineNextEvent() was called, 
        and that set the pair.activeDomain we use here. Two possibilities
        If active:
        * IV domain: 
            drawEventType returns either:
            - Reaction. Escape flag is not set (wouldn't be used anyway, 
              because only CoM is updated)
            - Escape. drawEventType sets escape flag in IV domain, fixes exit 
              point on IV shell, used via propagatePair.
        * CoM domain:
            drawEventType returns:
            - Escape. drawEventType sets escape flag in CoM domain, fixes exit 
              point on CoM shell, used via propagatePair.
	'''
	eventType = pair.activeDomain.drawEventType( pair.dt )
        log.info( '\tPair ' + str( eventType ) )


	if eventType == EventType.REACTION:
	    if len( pair.rt.products ) == 1:
		species3 = pair.rt.products[0]

		# calculate new R
		newCoM = pair.drawNewCoM( pair.dt )
		
		assert self.distance( pair.CoM, newCoM ) + species3.radius <\
		    pair.shellSize

		self.applyBoundary( newCoM )

		self.removeParticle( particle1 )
		self.removeParticle( particle2 )

		particle = self.createParticle( species3, newCoM )
		newsingle = self.createSingle( particle )

		self.reactionEvents += 1
		self.lastReaction = Reaction( pair.rt, [particle1, particle2],
					      [particle] )

                log.info( '\tNew %s. radius=%g. dt=%g.' % ( newsingle, newsingle.radius, newsingle.dt ) )

	    else:
		raise NotImplementedError,\
		      'num products >= 2 not supported.'

	    self.removeFromShellMatrix( pair )
	    pair.dt = -INF
	    return pair.dt

	# Escaping through a_r or escaping through a_R. Make use of escape flag 
	# magic.
	if eventType == EventType.ESCAPE:
	    self.propagatePair( pair )
	else:
	    raise SystemError, 'Bug: invalid eventType.'

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
	single1.pos, single2.pos = pair.drawNewPositions( self.t - pair.lastTime )
	self.applyBoundary(single1.pos)
	self.applyBoundary(single2.pos)
	assert self.checkOverlap( single1.pos, single1.getMinRadius(),
				  ignore = [ single1.particle, single2.particle] )
	assert self.checkOverlap( single2.pos, pair.single2.getMinRadius(),
				  ignore = [ single1.particle, single2.particle ] )
	assert pair.checkNewpos( single1.pos, single2.pos )
	self.moveParticle( single1.particle, single1.pos )
	self.moveParticle( single2.particle, single2.pos )
	single1.initialize( self.t )
	single2.initialize( self.t )
	self.addToShellMatrix( single1 ) # Add, not update.
	self.addToShellMatrix( single2 )
        log.info( '\tNew %s. radius=%g. dt=%g.' % ( single1, single1.radius, single1.dt ) )
        log.info( '\tNew %s. radius=%g. dt=%g.' % ( single2, single2.radius, single2.dt ) )


    def formInteractionOrPairOrMulti( self, single, neighbors ):
	assert neighbors
	bursted = []

        # Todo.
        # 1. Would be better to return surfaces seperately from getNeighbors.
        # 2. If multi within neighbors, return multi.
        # 3. Elif surface in surfacelist, try and return interaction.
        # 4. Else or if no interaction possible and single in surfacelist, try 
        # and return pair.
        # 5. Else. Multi.

        # Try interaction
	if isinstance( neighbors[0], Surface ):
	    obj = self.formInteraction( single, neighbors[0], neighbors[1:] )
	    if obj:
		return obj, neighbors[1:]
        elif isinstance( neighbors[0], Single ):
            # Try forming a Pair only if singles are on same surface.
            if single.surface == neighbors[0].surface:
                obj = self.formPair( single, neighbors[0], neighbors[1:] )
                if obj:
                    return obj, neighbors[1:]

	# Then, a Multi.
        neighborsString = ''
        for n in neighbors:
            neighborsString += str(n) + ', '
	log.debug( '\t\tDebug. Try to form Multi:%s +\n\t\t\t[ %s ].' % (single, neighborsString) )
	minShell = single.getMinRadius() * ( 1.0 + self.MULTI_SHELL_FACTOR )

        # Todo. Add surfaces to multi somehow.
        #neighbors = [ n for n in neighbors if not isinstance( n, Surface ) ]
	#if not neighbors:
        #    return None, bursted

	neighborDists = self.objDistanceArray( single.pos, neighbors )
	neighbors = [ neighbors[i] for i in 
		      ( neighborDists <= minShell ).nonzero()[0] ]

	if not neighbors:
            log.debug( '\t\tDebug. Multi not needed.' )
	    return None, bursted

	closest = neighbors[0]
	if isinstance( closest, Single ) or isinstance( closest, Surface ):

	    multi = self.createMulti()
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
	    return multi, bursted

	elif isinstance( closest, Multi ):

	    multi = closest
	    log.info( '\tmulti merge %s %s' % ( single, multi ) )

	    self.removeFromShellMatrix( multi )

	    self.addToMulti( single, multi )
	    self.removeFromShellMatrix( single )
	    for neighbor in neighbors[1:]:
		self.addToMultiRecursive( neighbor, multi )

	    multi.initialize( self.t )

            # Add 1 shell for each particle to main simulator's shellMatrix.
	    self.addToShellMatrix( multi )
	    self.updateEvent( self.t + multi.dt, multi )
	    return multi, bursted

	assert False, 'do not reach here'


    '''
    Find largest possible cylinder around particle and surface, such that it 
    is not interfering with other shells. Miedema's algorithm.
    '''
    def formInteraction( self, single, surface, bursted ):
	log.debug( '\t\tDebug. Try formInteraction: %s + %s.' % (single, surface) )

        particle = single.particle
	projectedPoint, projectionDistance = surface.projectedPoint( single.pos )

        # For interaction with a planar surface, decide orientation.  
	orientationVector = cmp(projectionDistance, 0) * surface.unitZ 
        particleDistance = abs( projectionDistance )

        particleRadius = particle.species.radius

        '''
        Initialize dr, dzl, dzr.

        For an interaction with a PlanarSurface:
        * dr is the radius of the cylinder.
        * dzl determines how much the cylinder is sticking out on the 
        other side of the surface, measured from the projected point.
        * dzr is the distance between the particle and the edge of the 
        cylinder in the z direction.

        For interaction with a CylindricalSurface:
        * dr is the distance between the particle and the edge of the 
        cylinder in the r direction.
        * dzl is the distance from the projected point to the left edge of 
          the cylinder.
        * dzr is the distance from the projected point to the right edge 
          of the cylinder.
        '''
        mindr  = particleRadius * ( 1.0 + self.SINGLE_SHELL_FACTOR )
        mindzl = particleRadius * ( 1.0 + self.SINGLE_SHELL_FACTOR )
        mindzr = particleRadius * ( 1.0 + self.SINGLE_SHELL_FACTOR )
        dr = self.getMaxShellSize()
	dzl = self.getMaxShellSize()
        dzr = self.getMaxShellSize()
	if isinstance( surface, Box ):
            dzr -= particleDistance
	elif isinstance( surface, Cylinder ):
            dr -= particleDistance

	allNeighbors = self.getNeighborsWithinRadiusNoSort( projectedPoint, dr, ignore=[single,] )

	for object in allNeighbors:
	    shell = object.shellList[0]
            # Also works if shell is a cylinder parallel to the same 
            # CylindricalSurface, because getRadius magically returns the size 
            # of the cylinder then.
	    objectRadius = shell.radius

            '''
            Make bursted singles look bigger, like in formPair, because the 
            size of their shell is only particle.radius (not yet multiplied by 
            SINGLE_SHELL_FACTOR) (and no we can not do that immediately after 
            they are bursted, singles might start overlapping).
            '''
            if object.dt == 0.0 and object.getD() > 0:
                # This is one of the bursted singles.
                # Or a particle that just escaped it's multi!!!
                # Should account for this also in formPair???
                objectRadius *= ( 1.0 + self.SINGLE_SHELL_FACTOR )
                #assert bursted.__contains__( object ), 'bursted=%s does not contain %s. radius=%g.' %(bursted, object, object.radius)

	    objectVector = shell.origin - projectedPoint

            # Calculate zi and ri for this object.
	    zi = numpy.dot( objectVector, orientationVector )
	    vectorZ = zi*numpy.array(orientationVector)
	    vectorR = numpy.array(objectVector) - numpy.array(vectorZ)
	    ri = numpy.linalg.norm( vectorR )

            # Calculate dri for this object.
	    dri = ri - objectRadius
	    if isinstance( surface, Cylinder ):
		dri -= particleDistance

            # Calculate dzli or dzri (both are usually positive values).
            if zi < 0:
                # Calculate dzli for this object.
                dzli = - zi - objectRadius

                # Miedema's algorithm left side.
                if dzli < dzl and dri < dr:
                    if dzli > dri:
                        dzli = dzli
                    else:
                        dr = dri
            else:
                # Calculate dzri for this object.
                dzri = zi - objectRadius

                if isinstance( surface, Box ):
                    # On the particle side (right side), do Miedema's 
                    # algorithm relative to the particle's position in the z 
                    # direction. 
                    dzri -= particleDistance

                # Miedema's algorithm right side.
                if dzri < dzr and dri < dr:
                    if dzri > dri:
                        dzi = dzri
                    else:
                        dr = dri

        if dr < mindr or dzl < mindzl or dzr < mindzr:
            log.debug( '\t\tDebug. Interaction not possible: %s + %s. dr=%g. mindr=%g. dzl=%g. mindzl=%g. dzr=%g. mindzr=%g.' % ( single, surface, dr, mindr, dzl, mindzl, dzr, mindzr ) )
            return None

        '''
        Compute radius and size of new cylinder.

        Todo. Should we make it even smaller? Like updateSingle: depending on 
        diffusion constant etc.
        '''
	if isinstance( surface, Box ):
            # On the other side (left side) than the particle's side, make 
            # sure there is just enough space for the particle to stick out, 
            # but not more.
            dzl = particleRadius 
	    radius = dr
            # sizeOfDomain is really different from size of cylinder!
            sizeOfDomain = particleDistance + dzr - surface.Lz
	    size = ( particleDistance + dzl + dzr ) / 2
	elif isinstance( surface, Cylinder ):
	    radius = dr + particleDistance
            # sizeOfDomain is size of cylinder * 2.
            sizeOfDomain = dzl + dzr
	    size = ( sizeOfDomain ) / 2

        # Compute new origin of cylinder.
        shiftZ = size - dzl
        origin = projectedPoint + shiftZ * orientationVector

        if isinstance( surface, Box ):
            # Compute new particle offset relative to origin of the domain in 
            # the z-direction (z=0), which is *at* the surface boundary. z=L 
            # is at the end of the cylinder.
            # (minRadius correction in the constructor).
            particleOffset = [ 0, particleDistance - surface.Lz ]
	elif isinstance( surface, Cylinder ):
            # Compute new particle offset in r and z-direction relative to 
            # origin of new cylinder. In the z-direction, the origin (z=0) is 
            # at the left boundary. z=L is at the right boundary.
            # (minRadius correction in the constructor).
            particleOffset = [ particleDistance, dzl ]


        '''
        Create interaction.
        '''
        assert single.dt == 0.0 and single.getMobilityRadius() == 0.0

        reactionTypes = self.getReactionType1( particle.species )
	interactionType = self.getInteractionType( particle.species, surface )


        interaction = surface.defaultInteractionSingle( particle, surface, reactionTypes, interactionType, origin, radius, orientationVector, size, particleOffset, projectedPoint, sizeOfDomain )
        interaction.initialize( self.t )

	self.addToShellMatrix( interaction )
	interaction.dt, interaction.eventType, interaction.activeDomain = interaction.determineNextEvent( )
	self.addSingleEvent( interaction )

	self.removeFromShellMatrix( single )

	cylinder = interaction.shellList[0]
	log.info( '\tNew %s. radius=%g. size=%g. dt=%g.' % (interaction, cylinder.radius, cylinder.size, interaction.dt) )

        assert self.checkObj( interaction )

	return interaction


    '''
    Decide if pair makes sense, and create it if so.
    '''
    def formPair( self, single1, single2, bursted ):
	log.debug( '\t\tDebug. Try formPair %s.' %
	           'Pair( %s, %s )' % ( single1.particle, 
	                                single2.particle ) )
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
	assert distanceFromSigma >= 0, '(pair gap) between %s and %s = %g < 0' \
	    % ( single1, single2, distanceFromSigma )

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
			    distanceFromSigma * 100 + sigma + shellSizeMargin )

	if minShellSizeWithMargin >= maxShellSize:
	    log.debug( '\t\tDebug. %s not formed: minShellSize >= maxShellSize.' %
		       ( 'Pair( %s, %s )' % ( single1.particle, 
					      single2.particle ) ) )
	    return None

	'''
	Here, we have to take into account of the bursted Singles in
	this step.  The simple check for closest below could miss
	some of them, because sizes of these Singles for this
	distance check has to include SINGLE_SHELL_FACTOR, while
	these bursted objects have zero mobility radii.  This is not
	beautiful, a cleaner framework may be possible.
	'''

	closest, closestShellDistance = DummySingle(), INF
	for b in bursted:
	    if isinstance( b, Single ):
		d = self.distance( com, b.pos ) \
		    - b.getMinRadius() * ( 1.0 + self.SINGLE_SHELL_FACTOR )
		if d < closestShellDistance:
		    closest, closestShellDistance = b, d

	if closestShellDistance <= minShellSizeWithMargin:
	    log.debug( '\t\tDebug. %s not formed: squeezed by bursted neighbor\n\t\t\t%s.' %
		       ( 'Pair( %s, %s )' % ( single1.particle, 
					      single2.particle ), closest ) )
	    return None

	c, d = self.getClosestObj( com, ignore=[ single1, single2 ] )
	if d < closestShellDistance:
	    closest, closestShellDistance = c, d

	log.debug( '\t\tDebug. Try pair. closest neighbor=%s.\n\t\t\tclosestShellDistance=%g. minShellWithMargin=%g.' %
		   ( closest, closestShellDistance, minShellSizeWithMargin ) )

	if isinstance( closest, Single ):

	    D_closest = closest.particle.species.D
	    D_tot = D_closest + D12
	    closestDistance = self.distance( com, closest.pos )

	    closestMinSize = closest.getMinRadius()
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
	    assert isinstance( closest, ( Pair, Multi, Surface, DummySingle ) ), closest

	    shellSize = closestShellDistance / SAFETY

	if shellSize <= minShellSizeWithMargin:
	    log.debug( '\t\tDebug. %s not formed: squeezed by %s.\n\t\t\tshellSize=%g. minShellSizeWithMargin=%g.' %
		       ( 'Pair( %s, %s )' % ( single1.particle, 
					      single2.particle ), closest, shellSize, minShellSizeWithMargin ) )
	    return None

	d1 = self.distance( com, single1.pos )
	d2 = self.distance( com, single2.pos )

	if shellSize < max( d1 + single1.getMinRadius() *
			    ( 1.0 + self.SINGLE_SHELL_FACTOR ), \
				d2 + single2.getMinRadius() * \
				( 1.0 + self.SINGLE_SHELL_FACTOR ) ) * 1.3:
	    log.debug( '\t\tDebug. %s not formed: singles are better.' %
		       'Pair( %s, %s )' % ( single1.particle, 
					    single2.particle ) )
	    return None

	# 3. Ok, Pair makes sense.  Create one.
	shellSize = min( shellSize, maxShellSize )

	pair = self.createPair( single1, single2, shellSize )

	self.removeFromShellMatrix( single1 )
	self.removeFromShellMatrix( single2 )
	self.addToShellMatrix( pair )


	# Formerly known as the impure function Pair.determineNextEvent().
	dtSingleReaction, reactingSingle = pair.drawSingleReactionTime( )
	dtEscape, activeDomain = pair.drawEscapeOrReactionTime( )

	if dtSingleReaction < dtEscape:
	    pair.eventType = EventType.REACTION # This is single (!) reaction.
	    pair.reactingSingle = reactingSingle
	    pair.dt = dtSingleReaction
	else:
	    pair.eventType = EventType.ESCAPE
	    pair.activeDomain = activeDomain
	    pair.dt = dtEscape

	pair.lastTime = self.t
	assert pair.dt >= 0

	self.addPairEvent( pair )
	# After returning, single1 is scheduled to past in fireSingle.
	self.removeEvent( single2 )

	assert closestShellDistance == INF or pair.shellSize < closestShellDistance
	assert pair.shellSize >= minShellSizeWithMargin
	assert pair.shellSize <= maxShellSize

	log.info( '\tNew %s. radius=%g. pairDistance=%g. dt=%g.' % ( pair, pair.shellSize, pairDistance, pair.dt ) )

	assert self.checkObj( pair )
	return pair


    ##########################################################################
    '''
    Methods for Multis.
    '''
    def createMulti( self ):
	multi = Multi( self )
	#multi.initialize( self.t )
	return multi


    def fireMulti( self, multi ):
	sim = multi.sim
	sim.step()
	#sim.sync()

	if sim.lastReaction:
	    log.info( '\tbd reaction' )

	    self.breakUpMulti( multi )
	    self.reactionEvents += 1
	    self.lastReaction = sim.lastReaction
	    return -INF

	if sim.escaped:
	    log.info( '\tmulti particle escaped.' )

	    self.breakUpMulti( multi )
	    return -INF

	#log.info( '\tmulti stepped %d steps, duration %g, dt = %g' %
	#          ( additionalSteps + 1, sim.t - startT + sim.dt, dt ) )
	return multi.dt


    def breakUpMulti( self, multi ):
	self.removeFromShellMatrix( multi )
	singles = []
	for particle in multi.sim.particleList:
	    single = self.createSingle( particle )
	    singles.append( single )
	return singles


    def burstMulti( self, multi ):
	#multi.sim.sync()
	singles = self.breakUpMulti( multi )
	return singles


    '''
    Add 'obj' to multi using addToMulti().
    
    If 'obj' is a Single, also add any neighbors that lie within 'shellSize' 
    of 'obj', where shellSize is MULTI_SHELL_FACTOR times the radius of the 
    particle.

    If 'obj' is itself a Multi, merge them.
    '''
    def addToMultiRecursive( self, obj, multi ):
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
	    shellSize = obj.particle.species.radius *\
		( 1.0 + self.MULTI_SHELL_FACTOR )
	    neighbors = self.getNeighborsWithinRadiusNoSort( obj.pos, shellSize,
							     ignore=[obj,] )
	    bursted = self.burstNonMultis( neighbors )
	    neighborDists = self.objDistanceArray( obj.pos, bursted )
	    neighbors = [ bursted[i] for i in 
			  ( neighborDists <= shellSize ).nonzero()[0] ]

	    for obj in neighbors:
		self.addToMultiRecursive( obj, multi )

	elif isinstance( obj, Multi ):
	    if not obj.sim.particleList[0] in multi.sim.particleList:
		self.mergeMultis( obj, multi )
		self.removeFromShellMatrix( obj )
		self.removeEvent( obj )
	    else:
		log.debug( '\t\tDebug. %s already added. skipping.' % obj )
	else:
	    assert False, 'do not reach here.'  # Pairs are bursted



    '''
    Add to multi:
        - the particle
        - a shell a bit bigger than radius (using MULTI_SHELL_FACTOR)
    '''   
    def addToMulti( self, single, multi ):
	log.info( '\tAdding %s to\n\t\t%s' % ( single, multi ) )

	shellSize = single.particle.species.radius * \
	    ( 1.0 + self.MULTI_SHELL_FACTOR )
	multi.addParticle( single.particle )
	multi.addShell( single.pos, shellSize )


    '''
    merge multi1 into multi2
    '''
    def mergeMultis( self, multi1, multi2 ):
	log.info( '\tmerging %s to %s' % ( multi1, multi2 ) )
	assert not multi1.sim.particleList[0] in multi2.sim.particleList

	for i, particle in enumerate( multi1.sim.particleList ):
	    # FIXME: shells should be renewed
	    multi2.addParticle( particle )
	    shell = multi1.shellList[i]
	    multi2.addShell( shell.origin, shell.radius )

	multi2.initialize( self.t )



