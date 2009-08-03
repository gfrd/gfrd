#!/usr/env python
import weakref
import math
import numpy

from gfrdbase import *

from single import *
from pair import *
from multi import *
from shape import Sphere, Cylinder

from itertools import izip
from vtklogger import VTKLogger

class Delegate( object ):

    def __init__( self, obj, method ):
        self.obj = weakref.proxy( obj )
        self.method = method

    def __call__( self, *arg ):
        return self.method( self.obj, *arg )


# Notes:
# EGFRDSimulator.distance takes into account the periodic boundary conditions.  
# Very handy!
# Returning a dt to the scheduler reschedules the event, -INF removes it.
class EGFRDSimulator( ParticleSimulatorBase ):
    def __init__( self, logdir = None ):
        if logdir:
            log.info('Started vtk logger: ' + logdir)
            self.vtklogger = VTKLogger( self, logdir )
        else:
            self.vtklogger = None

        self.sphereMatrix = SphereMatrix()
        self.cylinderMatrix = CylinderMatrix()
        self.boxMatrix = BoxMatrix()
        self.objectMatrices = [ self.sphereMatrix, self.cylinderMatrix, self.boxMatrix ]
        #self.sm2 = pObjectMatrix()

        ParticleSimulatorBase.__init__( self )

        self.MULTI_SHELL_FACTOR = 0.05
        self.SINGLE_SHELL_FACTOR = 1.0  # 0.1

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
        self.lastEvent = None
        self.lastReaction = None
        self.isDirty = True
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
        for p in particles:
            # Get reference to real Particle. I think this is needed.
            particle = Particle( p.species, serial=p.serial, surface=p.surface )
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


    def stop( self, t ):
        if self.vtklogger:
            self.vtklogger.stop()
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
    These methods returns a tuple ( neighbors, distances ).
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
        distanceToSurface, closestSurface = self.getClosestSurface( pos ) 
        if distanceToSurface < closestDistance:
            closestDistance = distanceToSurface
            closestObject = closestSurface

        return closestObject, closestDistance


    '''
    Returns sorted list of pairs:
    - distance to surface
    - surface itself

    We can not use objectmatrix, it would miss a surface if the origin of the 
    surface would not be in the same or neighboring cells as pos.
    '''
    def getClosestSurface( self, pos ):
        surfaces = [ None ]
        distances = [ INF ]
        for surface in self.surfaceList:
            posTransposed = cyclicTranspose( pos, surface.outside.origin, self.worldSize )
            distances.append( surface.signedDistanceTo( posTransposed ) )
        return min( zip( distances, surfaces ))


    '''
    Todo: no surface detection needed?
    '''
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
        distanceToSurface, closestSurface = self.getClosestSurface( pos ) 
        if distanceToSurface < radius:
            neighbors.append( closestSurface )
            distances.append( distanceToSurface )
        elif distanceToSurface < closestDistance:
            closestDistance = distanceToSurface
            closestSingle = closestSurface

        return neighbors, distances, closestSingle, closestDistance


    # Todo: do better for cylinders, using periodic boundary condition.
    def objDistance( self, pos, obj ):
        dists = numpy.zeros( len( obj.shellList ) )
        for i, shell in enumerate( obj.shellList ):
            try:
                shell.orientation
                size = math.sqrt( pow(shell.radius, 2) + pow(shell.size, 2))
            except:
                size = shell.radius
                
            dists[i] = self.distance( pos, shell.origin ) - size
        return min( dists )


    def objDistanceArray( self, pos, objs ):
        dists = numpy.array( [ self.objDistance( pos, obj ) for obj in objs ] )
        """
        log.debug("Todo: better periodic boundary condition handling.")
        dists = numpy.array( [ obj.distanceTo( pos ) for obj in objs ] )
        """
        return dists
            

    ##########################################################################
    '''
    consistency checkers
    Todo.
    '''
    def checkObj( self, obj ):
        obj.check()

        allshells = [ ( obj, i ) for i in range( len( obj.shellList ) ) ]
        for i, shell in enumerate( obj.shellList ):
            closest, distance = self.getClosestObj( shell.origin,
                                                    ignore = [obj] )
            radius = shell.radius
            assert radius <= self.getUserMaxShellSize(),\
                '%s shell radius larger than user-set max shell radius' % \
                str( ( obj, i ) )
            assert radius <= self.getMaxShellSize(),\
                '%s shell radius larger than simulator cell radius / 2' % \
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
	single = particle.surface.defaultSingle( particle, rt, Delegate( self, EGFRDSimulator.distance ) )
	single.initialize( self.t )
	self.addToShellMatrix( single )
	self.addSingleEvent( single )
	return single


    def fireSingle( self, single ):
	assert abs( single.dt + single.lastTime - self.t ) <= 1e-18 * self.t

	# Reaction.
	if single.eventType == EventType.REACTION:
	    log.info( 'single reaction %s' % str( single ) )

	    self.propagateSingle( single )
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

	'''
	Propagate, if not reaction.
	Handle immobile case first.
	'''
	if single.getD() == 0:
	    # no propagation, just calculate next reaction time.
	    # (!) Not a pure function.
	    single.dt, single.eventType, single.activeDomain = single.determineNextEvent( ) 
	    single.lastTime = self.t
	    return single.dt
	
	'''
	Propagate this particle to the exit point on the shell. This is done by 
	setting the escape flag in the active domain.
	This is the only time we need this activeDomain stuff.
	'''
	single.activeDomain.drawEventType( single.dt )
	self.propagateSingle( single )

	# (2) Clear volume.
	minShell = single.getMinRadius() * ( 1.0 + self.SINGLE_SHELL_FACTOR )
	closeNeighbors, distances, closest, closestShellDistance = \
		self.getNeighbors( single.pos, minShell, ignore=[single,] )

	bursted = []
	if closeNeighbors:
	    bursted = self.burstNonMultis( closeNeighbors )
	    obj, b = self.formInteractionOrPairOrMulti( single, bursted )
	    bursted.extend( b )
	    if obj:
		single.dt = -INF # remove by rescheduling to past.
		return single.dt

	    # if nothing was formed, recheck closest and restore shells.
	    closest, closestShellDistance = \
		self.getClosestObj( single.pos, ignore = [ single, ] )

	# All neighbors are more than minShell away.
	self.updateSingle( single, closest, closestShellDistance )
	bursted = uniq( bursted )
	burstedSingles = [ s for s in bursted if isinstance( s, Single ) ]
	self.restoreSingleShells( burstedSingles )
	log.info( 'single shell pos=(%3.1f %3.1f %3.1f) radius=%g dt=%g.' % ( single.pos[0], single.pos[1], single.pos[2], single.radius, single.dt ) )
	return single.dt


    def burstSingle( self, single ):
	assert self.t >= single.lastTime
	assert self.t <= single.lastTime + single.dt
	assert single.radius >= single.getMinRadius()

	oldpos = single.pos
	mobilityRadius = single.getMobilityRadius()
	self.propagateSingle( single )
	assert self.distance( single.pos, oldpos ) <= mobilityRadius

	self.updateEvent( self.t, single )


    ''' 
    The difference between a burst and a propagate is that a burst always takes 
    place before the actual scheduled event for the single, while propagateSingle 
    can be called for an escape event.
    Another subtle difference is that burstSingle always reschedules (updateEvent) 
    the single, while just calling propagateSingle does not. This works if 
    single.dt is returned to the scheduler after calling propagateSingle.
    '''
    def propagateSingle( self, single ):
	single.pos = single.drawNewPosition(self.t - single.lastTime) 
	self.applyBoundary( single.pos )
	assert self.checkOverlap( single.pos, single.getMinRadius(),
				  ignore = [ single.particle ] )
	self.updateOnParticleMatrix( single.particle, single.pos )
	single.initialize( self.t )
	self.updateShellMatrix( single )


    def restoreSingleShells( self, singles ):
	for single in singles:
	    assert single.isReset()
	    c, d = self.getClosestObj( single.pos, ignore = [single,] )

	    self.updateSingle( single, c, d )
	    self.updateEvent( self.t + single.dt, single )
	    log.debug( 'restore shell %s %g dt %g closest %s %g' %
		       ( single, single.radius, single.dt, c, d ) )


    # Draw new shell + new event time.
    def updateSingle( self, single, closest, distanceToShell ): 
	if isinstance( closest, Single ): # and isinstance( closest.shellList[0], Sphere)
	    distanceToClosest = self.distance( single.pos, closest.pos )
	    shellSize = self.calculateSingleShellSize( single, closest, 
						       distanceToClosest,
						       distanceToShell )
	else:  # Pair or Multi
	    shellSize = distanceToShell / SAFETY
	    shellSize = max( shellSize, single.getMinRadius() )

	shellSize = min( shellSize, self.getMaxShellSize() )
	single.radius = shellSize

	single.dt, single.eventType, single.activeDomain = single.determineNextEvent( )
	single.lastTime = self.t
	# No need for self.updateEvent(), single.dt is returned from fireSingle().
	self.updateShellMatrix( single )


    def calculateSingleShellSize( self, single, closest, 
				  distance, shellDistance ):
	minSize1 = single.getMinRadius()
	D1 = single.getD()
	if D1 == 0:
	    return minSize1

	D2 = closest.getD()
	minSize2 = closest.getMinRadius()
	minSize12 = minSize1 + minSize2
	sqrtD1 = math.sqrt( D1 )
	shellSize = min( sqrtD1 / ( sqrtD1 + math.sqrt( D2 ) )
			 * ( distance - minSize12 ) + minSize1,
			 shellDistance )
	shellSize /= SAFETY
	shellSize = max( shellSize, minSize1 ) # not smaller than the radius
	return shellSize


    def fireSingleReaction( self, single, interactionType = None ):
	reactantSpecies = single.particle.species
	oldpos = single.particle.pos.copy()
	
	rt = single.drawReactionType()
	
	if isinstance( rt, SurfaceUnbindingReactionType ):
	    raise NotImplementedError, 'SurfaceUnbindingReactionType!'
	if interactionType:
	    pass
	    #'nu in fireSingleReaction voor cylinder'
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

	    currentSurface = single.particle.surface
	    self.removeParticle( single.particle )
	    newparticle = self.createParticle( productSpecies, oldpos, currentSurface )
	    newsingle = self.createSingle( newparticle )

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
		log.info( 'no space for product particles.' )
		raise NoSpace()

	    currentSurface = single.particle.surface
	    self.removeParticle( single.particle )
	    particle1 = self.createParticle( productSpecies1, newpos1, currentSurface )
	    particle2 = self.createParticle( productSpecies2, newpos2, currentSurface )
	    newsingle1 = self.createSingle( particle1 )
	    newsingle2 = self.createSingle( particle2 )

	    self.lastReaction = Reaction( rt, [single.particle], 
					  [particle1, particle2] )

	    log.info( 'products; %s %s' % 
		      ( str( newsingle1 ), str( newsingle2 ) ) )

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
	assert single1.particle.surface == single2.particle.surface
	rt = self.reactionTypeMap2.get( (single1.particle.species, 
	    single2.particle.species) )
	pair = single1.particle.surface.defaultPair(single1, single2, shellSize, 
		rt, Delegate( self, EGFRDSimulator.distance ), self.getWorldSize())
	pair.initialize( self.t )
	return pair


    def firePair( self, pair ):
	assert self.checkObj( pair )
	log.info( 'fire: %s eventType %s' % ( pair, pair.eventType ) )

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

	'''
	1. All other cases
	
	Decide if this is a *pair* reaction (0) or an escape (1) through either r 
	or R. 
	
	Reaction: IV domain is active, but escape flag is not set (wouldn't be 
	used anyway, because only CoM is updated)
	Escape: escape flag is set in active domain, fixes exit point on 
	appropriate shell.
	'''
	eventType = pair.activeDomain.drawEventType( pair.dt )

	if eventType == EventType.REACTION:
	    log.info( 'reaction' )

	    if len( pair.rt.products ) == 1:
		species3 = pair.rt.products[0]

		# calculate new R
		newCoM = pair.drawNewCoM( pair.dt )
		
		assert self.distance( pair.CoM, newCoM ) + species3.radius <\
		    pair.shellSize

		currentSurface = particle1.surface
		self.applyBoundary( newCoM )

		self.removeParticle( particle1 )
		self.removeParticle( particle2 )

		particle = self.createParticle( species3, newCoM, currentSurface  )
		newsingle = self.createSingle( particle )

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

	# Escaping through a_r or escaping through a_R. Make use of escape flag 
	# magic.
	if eventType == EventType.ESCAPE:
	    log.debug( 'pairDistance = %g, dt = %g, %s' %
			   ( pair.pairDistance, pair.dt, pair.pgf.dump() ) )
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
	self.updateOnParticleMatrix( single1.particle, single1.pos )
	self.updateOnParticleMatrix( single2.particle, single2.pos )
	single1.initialize( self.t )
	single2.initialize( self.t )
	self.addToShellMatrix( single1 ) # Add not update.
	self.addToShellMatrix( single2 )


    def formInteractionOrPairOrMulti( self, single, neighbors ):
	assert neighbors
	bursted = []

	# Try interaction
	if isinstance( neighbors[0], Surface ):
	    obj = self.formInteraction( single, neighbors[0], neighbors[1:] )
	    if obj:
		return obj, neighbors[1:]


	# Try forming a Pair.
	if isinstance( neighbors[0], Single ): # and isinstance( neighbors[0].shellList[0], Sphere)
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


    # Todo.
    # Find largest possible cylinder around particle, such that it is not
    # interfering with other particles. Miedema's algorithm.
    #def formCylinder( self, single, surface, interactionType ):
    def formInteraction( self, single, surface, bursted ):
	# Todo.
	assert bursted == []
	log.debug( 'formInteraction: %s' % (single) )

	orientation = surface.orientationZ
	origin, posVector, r0 = surface.calculateProjectionVectors( single.pos )

	# Todo. Applyboundary.
	assert numpy.linalg.norm(single.pos - origin) >= surface.radius

	dr = self.getMaxShellSize() #INF #CYLINDRICAL_SHELL_MAX_WIDTH - r0
	dz = self.getMaxShellSize() #INF #CYLINDRICAL_SHELL_MAX_HALF_LENGTH

	#cylinder = surface.interactionSingle( origin, 

	allNeighbors = getNeighborsWithinRadiusNoSort( origin, INF, ignore=[single,] )

	for object in allNeighbors:
	    rhoi = object.shellList[0].radius

	    objectVector = object.shellList[0].origin - origin

	    # Calculate dz for this object.
	    zi = numpy.dot( objectVector, orientation )
	    dzi = abs(zi) - rhoi
	    temp = zi*numpy.array(orientation)
	    temp2 = numpy.array(objectVector) - numpy.array(temp)

	    # Calculate dr for this object.
	    ri = numpy.linalg.norm( temp2 )
	    dri = ri - r0 - rhoi

	    if dzi < dz and dri < dr:
		if dzi > dri:
		    dz = dzi
		else:
		    dr = dri

	assert dr > 0
	bursted = uniq( bursted )
	burstedSingles = [ s for s in bursted if isinstance( s, Single ) ]
	# Probably everything would work just fine still if we didn't restore
	# here, since those burstedSingles are already in the scheduler with a
	# dt=0 because they have been given a shell with size minSize.
	self.restoreSingleShells( burstedSingles )

	radius_a = surface.radius
	radius_b = r0 + dr
	halfLength = abs(dz)
	# Make cylinder with radius b and half-length dz. 
	cylinder = self.createCylinder( single, origin, orientation, radius_a, r0, radius_b, halfLength, interactionType )

	self.removeFromShellMatrix( single )

	self.shellMatrix.addCylinder( (cylinder, 0), cylinder.shellList[0] )
	cylinder.determineNextEvent( self.t )
	self.addCylinderEvent( cylinder )
	log.info( 'cylinder shell pos %s radius %g halfLength %g dt %g' % (cylinder.origin, cylinder.radius, cylinder.halfLength, cylinder.dt) )
	#Todo:
	#assert self.checkObj( cylinder )
	return cylinder


    # Todo.
    def formInteraction( self, single, surface, bursted ):
	# Todo.
	assert bursted == []


	

	'''
	Here, we have to take into account of the bursted Singles in
	this step.  The simple check for closest below could miss
	some of them, because sizes of these Singles for this
	distance check has to include SINGLE_SHELL_FACTOR, while
	these bursted objects have zero mobility radii.  This is not
	beautiful, a cleaner framework may be possible.
	'''
	# TODO 
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


    # Decide if pair makes sense.
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
	    log.debug( '%s not formed: minShellSize >= maxShellSize' %
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

	log.info( '%s, dt=%g, pairDistance=%g, shell=%g,' %
		  ( pair, pair.dt, pairDistance, pair.shellSize ) + 
		  'closest=%s, closestShellDistance=%g' %
		  ( closest, closestShellDistance ) )

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
	    singles.append( single )
	return singles


    def burstMulti( self, multi ):
	#multi.sim.sync()
	singles = self.breakUpMulti( multi )
	return singles


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
	    multi2.addShell( shell.origin, shell.radius )

	multi2.initialize( self.t )



