
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
        self.addToShellMatrix( single )
        self.addSingleEvent( single )
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
        multi2.addShell( shell.origin, shell.size )

    multi2.initialize( self.t )


