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
        

    def initialize( self ):
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
            self.addToShellMatrix( single )
            self.addSingleEvent( single )

        nParticles = sum( [ s.pool.size for s in self.speciesList.values() ] )
        assert nParticles == len( particles )

        self.isDirty = False


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


    def step( self ):
        if self.vtklogger:
            self.vtklogger.log()

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


    def burstNonMultis( self, neighbors ):
        bursted = []

        for obj in neighbors:
            if not isinstance( obj, Multi ):
                b = self.burstObj( obj )
                bursted.extend( b )
            else:
                bursted.append( obj )

        return bursted


    def clearVolume( self, pos, radius, ignore=[] ):
        neighbors = self.getNeighborsWithinRadiusNoSort( pos, radius, ignore )

        return self.burstObjs( neighbors )



    # Todo: do better for cylinders, using periodic boundary condition.
    def objDistance( self, pos, obj ):
        dists = numpy.zeros( len( obj.shellList ) )
        for i, shell in enumerate( obj.shellList ):
            try:
                shell.orientation
                size = math.sqrt( pow(shell.radius, 2) + pow(shell.size, 2))
            except:
                size = shell.size
                
            dists[i] = self.distance( pos, shell.origin ) - size
        return min( dists )


    def objDistanceArray( self, pos, objs ):
        dists = numpy.array( [ self.objDistance( pos, obj ) for obj in objs ] )
        """
        log.debug("Todo: better periodic boundary condition handling.")
        dists = numpy.array( [ obj.distanceTo( pos ) for obj in objs ] )
        """
        return dists
            


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



    execfile("/home/miedema/code/epdp-0.3/egfrd-single.py")
    execfile("/home/miedema/code/epdp-0.3/egfrd-pair.py")
    execfile("/home/miedema/code/epdp-0.3/egfrd-multi.py")
    execfile("/home/miedema/code/epdp-0.3/egfrd-matrix.py")
    execfile("/home/miedema/code/epdp-0.3/egfrd-check.py")
