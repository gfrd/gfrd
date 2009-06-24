
def createPair( self, single1, single2 ):
    assert single1.dt == 0.0
    assert single2.dt == 0.0
    assert single1.getMobilitySize() == 0.0
    assert single2.getMobilitySize() == 0.0

    species1 = single1.particle.species
    species2 = single2.particle.species
    rt = self.reactionTypeMap2.get( ( species1, species2 ) )

    pair = Pair( single1, single2, rt, 
                 Delegate( self, EGFRDSimulator.distance ),
                 self.getWorldSize() )
    pair.initialize( self.t )

    return pair


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


