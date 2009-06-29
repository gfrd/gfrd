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
            if not isinstance( obj, Multi ):
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
        else:  # Multi
            bursted = self.burstMulti( obj )
            self.removeEvent( obj )
            return bursted


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


    def getMatrixCellSize( self ):
        assert self.sphereMatrix.cellSize == self.cylinderMatrix.cellSize
        return self.sphereMatrix.cellSize


    def addToShellMatrix( self, obj, update=False ):
        for i, shell in enumerate( obj.shellList ):
            key = (obj, i)
            if isinstance( shell, Sphere ):
                self.sphereMatrix.add( key, shell.origin, shell.radius, update )
            elif isinstance( shell, Cylinder ):
                self.cylinderMatrix.add( key, shell, update )
            else: raise KeyError, 'Objecttype does not exit'


    def removeFromShellMatrix( self, obj ):
        for i, shell in enumerate( obj.shellList ):
            key = (obj, i)
            '''
            Check for type of shell -> remove key (is that the way to do it?)
            '''
            if isinstance( shell, Sphere ):
                self.sphereMatrix.remove( key )
            elif isinstance( shell, Cylinder ):
                self.cylinderMatrix.remove( key )
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
        closestSingle = DummySingle()
        closestDistance = INF

        for objectMatrix in [ self.sphereMatrix, self.cylinderMatrix ]:
            keys, distances = objectMatrix.getNeighbors( pos )
            '''
            Don't be clever and think you can do:
                keys, distances = objectMatrix.getNeighbors( pos, 1 )
            Because then you are ignoring the ignore list.
            '''
            for i, key in enumerate( keys ):
                single = key[0]
                if single not in ignore and distances[i] < closestDistance:
                    closestSingle, closestDistance = single, distances[i]
                    '''
                    Found yet a closer single. Break out of inner for loop 
                    and check other objectMatrices.
                    '''
                    break   
        return closestSingle, closestDistance


    def getNeighborsWithinRadiusNoSort( self, pos, radius, ignore=[] ):
        keys, _ =\
            self.sphereMatrix.getNeighborsWithinRadiusNoSort( pos, radius )
        keys2, _ =\
            self.cylinderMatrix.getNeighborsWithinRadiusNoSort( pos, radius )

        neighbors = uniq( [ key[0] for key in keys if key[0] not in ignore ] )
        neighbors2 = uniq( [ key[0] for key in keys2 if key[0] not in ignore ] )

        neighbors.extend( neighbors2 )
        return neighbors


    '''
    Returns all singles and the distances towards their *shell* if that shell is 
    within a distance radius of pos, plus the closest single and the distance 
    towards it shell that is just outside of radius.
    '''
    def getNeighbors( self, pos, radius=INF, ignore=[] ):
        # Diccionaries are more efficient for lookups.
        seen = dict.fromkeys( ignore )
        neighbors = []
        distances = []

        closestSingle = DummySingle()
        closestDistance = INF

        for objectMatrix in [ self.sphereMatrix, self.cylinderMatrix ]:
            keys, dists = objectMatrix.getNeighbors( pos )

            for i, key in enumerate( keys ):
                single = key[0]
                if not single in seen:
                    '''
                    Since a single can in theory have more than 1 shell, and 
                    for each shell there is an entry in the objectMatrix, we 
                    signal here that we have already seen this single.
                    This is a bug: seen[ key ] = None
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



    execfile("/home/miedema/code/epdp-0.3/egfrd-single.py")
    execfile("/home/miedema/code/epdp-0.3/egfrd-pair.py")
    execfile("/home/miedema/code/epdp-0.3/egfrd-multi.py")

