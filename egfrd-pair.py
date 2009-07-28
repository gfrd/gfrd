
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

