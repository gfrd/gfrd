
def moveSingle( self, single, pos ):
    single.pos = pos
    self.updateOnParticleMatrix( single.particle, pos )


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
    #print 'updateSingle'
    self.updateShellMatrix( single )


def createSingle( self, particle ):
    # Poissonian reactions.
    rt = self.getReactionType1( particle.species )
    # The type of defaultSingle depends on the surface this particle is 
    # on. Usually SphericalSingle or CylindricalSingle.
    single = particle.surface.defaultSingle( particle, rt, Delegate( self, EGFRDSimulator.distance ) )
    single.initialize( self.t )
    return single


def burstSingle( self, single ):
    assert self.t >= single.lastTime
    assert self.t <= single.lastTime + single.dt
    assert single.size >= single.getMinSize()

    # Set all domains to inactive.
    for domain in single.domains:
        domain.active = False

    oldpos = single.pos
    newpos = single.burst(self.t - single.lastTime)
    self.applyBoundary( newpos )

    assert self.distance( newpos, oldpos ) <= single.getMobilitySize()
    assert self.checkOverlap( newpos, single.particle.species.radius,\
                              ignore = [ single.particle, ] )

    # Move particle.
    self.moveSingle( single, newpos )

    single.initialize( self.t )
    self.updateShellMatrix( single )

    # Todo: is this not not needed for propagate?
    self.updateEvent( self.t, single )


def propagateSingle( self, single ):
    assert abs( single.dt + single.lastTime - self.t ) <= 1e-18 * self.t
    newpos = single.propagate(self.t - single.lastTime) 
    self.applyBoundary( newpos )

    assert self.checkOverlap( newpos, single.getMinSize(),
                              ignore = [ single.particle, ] )

    # Move particle.
    self.moveSingle( single, newpos )

    single.initialize( self.t )
    self.updateShellMatrix( single )


def fireSingleReaction( self, single, interactionType = None ):
    reactantSpecies = single.particle.species
    oldpos = single.particle.pos.copy()
    
    rt = single.drawReactionType()
    
    if isinstance( rt, SurfaceUnbindingReactionType ):
        raise NotImplementedError, 'SurfaceUnbindingReactionType!'

    if interactionType:
        print 'nu in fireSingleReaction voor cylinder'

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

        currentSurface = single.particle.surface
        self.removeParticle( single.particle )

        particle1 = self.createParticle( productSpecies1, newpos1, currentSurface )
        particle2 = self.createParticle( productSpecies2, newpos2, currentSurface )
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


def fireSingle( self, single ):
    # Reaction.
    if single.eventType == EventType.REACTION:

        log.info( 'single reaction %s' % str( single ) )

        # Todo: does it matter that we call updateEvent from here, while 
        # this wasn't called when this said propagateSingle?
        self.burstSingle( single )

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
    
    self.propagateSingle( single )

    # (2) Clear volume.

    minShell = single.getMinSize() * ( 1.0 + self.SINGLE_SHELL_FACTOR )

    closeNeighbors, distances, closest, closestShellDistance = \
            self.getNeighbors( single.pos, minShell, ignore=[single,] )

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
        
    log.info( 'single shell pos=(%3.1f %3.1f %3.1f) size=%g dt=%g.' % ( single.pos[0], single.pos[1], single.pos[2], single.size, single.dt ) )

    return single.dt
    # Reaction
    if shell.eventType == EventType.REACTION:
        self.burstShell( shell )
        self.removeFromShellMatrix( single )
        self.fireSingleReaction( shell )

        shell.dt = -INF  # remove this Shell from the Scheduler
        return shell.dt

    # Escape
    else:
        self.propagateShell( shell ) 
        pass


def calculateSingleShellSize( self, single, closest, 
                              distance, shellDistance ):

    minSize1 = single.getMinSize()
    D1 = single.getD()

    if D1 == 0:
        return minSize1

    #print 'closest = ', closest
    #print 'distance = ', distance
    #print 'shellDistance = ', shellDistance
    #log.debug( 'closest.size = %g' % (closest.size * 1e9) )
    #log.debug( 'distanceToClosestParticle - closest.size = %g' % ((distance - closest.shellList[0].size) * 1e9) )
    # Not necessarily true: see test7 object_matrix_test.py.
    #if distance < shellDistance:
    #    if self.vtklogger:
    #        self.vtklogger.stop()
    #     raise Exception, "stop"
    #assert shellDistance == distance - closest.shellList[0].size
    D2 = closest.getD()
    minSize2 = closest.getMinSize()
    minSize12 = minSize1 + minSize2
    sqrtD1 = math.sqrt( D1 )
        
    shellSize = min( sqrtD1 / ( sqrtD1 + math.sqrt( D2 ) )
                     * ( distance - minSize12 ) + minSize1,
                     shellDistance )
    #log.debug( 'shellSize = %g' % shellSize )
    shellSize /= SAFETY
    #log.debug( 'shellSize after safety = %g' % shellSize )
    shellSize = max( shellSize, minSize1 ) # not smaller than the radius

    return shellSize


def restoreSingleShells( self, singles ):

    for single in singles:
        assert single.isReset()
        c, d = self.getClosestObj( single.pos, ignore = [single,] )

        self.updateSingle( single, c, d )
        self.updateEvent( self.t + single.dt, single )
        log.debug( 'restore shell %s %g dt %g closest %s %g' %
                   ( single, single.size, single.dt, c, d ) )


