
# consistency checkers
# Todo.
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

