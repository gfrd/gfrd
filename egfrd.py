from egfrd1 import *
from itertools import izip

class EGFRDSimulator( EGFRDSimulator1 ):

    def __init__( self ):
        EGFRDSimulator1.__init__( self )


    def createSingle( self, particle ):
        # Poissonian reactions.
        rt = self.getReactionType1( particle.species )
        # The type of defaultSingle depends on the surface this particle is 
        # on. Usually SphericalSingle or CylindricalSingle.
        single = particle.surface.defaultSingle( particle, rt, Delegate( self, EGFRDSimulator.distance ) )
        single.initialize( self.t )
        return single


    # Returns sorted list with elements:
    # - distance to surface
    # - surface itself
    # - interactionType (only 1) this subSpace can have with that surface.
    def getClosestSurface( self, pos ):
        surfaces = ( s for s in self.surfaceList )
        #interactions = ( s[1] for s in self.surfaceList )
        try:
            # Todo: is this distance boundaries-safe?
            distances = ( (s.distance( pos ) ) for s in self.surfaceList )
            return min( izip( distances, surfaces ))
        except ValueError:
            return INF, None


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

        closeNeighbors, distances = self.getNeighbors( single.pos, minShell,
                                                       ignore=[single,] )

        # This is a bit tricky, but the last one in closeNeighbors
        # is the closest object to this Single.
        # getNeighbors() returns closeNeighbors within minShell *plus* one.
        closest = closeNeighbors.pop()
        closestShellDistance = distances[-1]

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
            
        log.info( 'single shell %s %g dt %g.' % ( single.pos, single.size, single.dt ) )

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

        assert self.checkOverlap( newpos, single.getMinSize(), \
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


# GRAVEYARD / Todo.
    """
    # fireSingle is more or less surface proof.
    ### A single particle either reacts or leaves it's shell. Here these
    ### events are executed.
    def fireSingle( self, single ):
        log.debug( '    %s.fireSingle: %s' % ('egfrd2', single) )

        ### CASE Reaction
        # Reaction.
        if single.eventType == EventType.REACTION:

            log.info( 'single reaction %s' % str( single ) )
            r = single.drawR( single.dt )

            self.propagateSingle( single )

            try:
                ### Remove the shell here. The actual particle will be removed
                ### in fireSingleReaction().
                self.removeFromShellMatrix( single )
                self.fireSingleReaction( single )
            except NoSpace:
                log.info( 'single reaction; placing product failed.' )
                self.addToShellMatrix( single )
                self.mainSimulator.rejectedMoves += 1
                single.reset()
                return single.dt

            ### Todo: how does this work?
            single.dt = -INF  # remove this Single from the Scheduler
            ### Particle removed in fireSingleReaction.
            return single.dt

        ### CASE Escape
        # Propagate, if not reaction.

        # Handle immobile case first.
        if single.getD() == 0:
            # no propagation, just calculate next reaction time.
            single.determineNextEvent( self.mainSimulator.t ) 
            return single.dt
        
        # Propagate this particle to the exit point on the shell.
        self.propagateSingle( single )

        # (2) Clear volume.

        ### The single was just propagated and initialized, so it's shell has
        ### size getMinSize().
        minShell = single.getMinSize() * ( 1.0 + self.SINGLE_SHELL_FACTOR )

        # We already know there are no other objects within getMinSize()
        # (particle radius), because propagateSingle checked it. Now check if
        # there are any within minShell, and also get closest outside of
        # minShell.
        closeNeighbors, distances = self.getNeighbors( single.pos, minShell,
                                                       ignore=[single,] )

        # This is a bit tricky, but the last one in closeNeighbors
        # is the closest object outside of this Single's shell.
        # getNeighbors() returns closeNeighbors within minShell *plus* one.
        # That's why we can not use getNeighborsWithinRadius.
        closest = closeNeighbors.pop()
        closestShellDistance = distances[-1]
        # So closest is now always the first object outside of minShell. If
        # closeNeigbors still contains objects, those lie inside of minShell,
        # and we should make a pair or multi with them.


        # Check surfaces.
        distanceToSurface, closestSurface = self.getClosestSurface( single.pos )
        # Possibilities, objects within minShell:
        # 1. object(s) + surface: try a pair, fail, try a multi (including the
        # surface).
        # 2. object(s) only: try a pair, if it fails, try a multi.
        # 3. surface only: try a shell that includes part of the surface.
        # 4. nothing: single (either to closest or to surface, whichever is
        # nearer).

        bursted = []
        
        if closeNeighbors:
            ### If a closeNeighbor is already a multi, don't burst it, but
            ### let it absorb the single in the next step.
            bursted = self.burstNonMultis( closeNeighbors )
            if distanceToSurface < minShell:
                log.debug( '    %s.fireSingle Case 1: neighbors + surface nearby' % ('egfrd2' ) )
                # Case 1.
                # Todo. Make this work in formPairOrMulti.
                # Todo. InteractionType.
                bursted.append(closestSurface)
            else:
                log.debug( '    %s.fireSingle Case 2: neighbors nearby, no surface' % ('egfrd2' ) )
                # Case 2.
                pass
            obj, b = self.formPairOrMulti( single, bursted )
            # obj can be a pair or a multi, b contains all singles that were
            # not added to it.
            # Why are we extending this list? It's already complete isn't it?
            bursted.extend( b )

            if obj:
                # Why are we returning here? formPairOrMulti returns all other
                # neighbors after a pair is formed with the first neighbor.
                # Shouldn't we restore those, I think we really should.
                # Same yields for bursted singles not added to the multi.
                #
                # Probably the answer to this question is that bursted singles
                # are give a size of minSize. So they are placed in the
                # scheduler and will be given a new shell soon. We could do it
                # here, but it isn't necessary.
                single.dt = -INF # remove by rescheduling to past.
                return single.dt

            # if nothing was formed, recheck closest and restore shells.
            # Correct. Maybe some other particle has come closer during the
            # bursting (while no multi had to be formed).
            # But again, we should also do this if a pair was formed.
            closest, closestShellDistance = \
                self.getClosestObj( single.pos, ignore = [ single, ] )
        else:
            if distanceToSurface < minShell:
                # Case 3.
                print 'Cylinder!'
                log.debug( '    %s.fireSingle Case 3: surface nearby, no neighbors' % ('egfrd2' ) )
                if isinstance( closestSurface, CylindricalSurface  ):
                    cylinder = self.formCylinder( single, closestSurface, True )
                    single.dt = -INF # remove by rescheduling to past.
                    # Return single.dt, not cylinder.dt!
                    return single.dt
                else:
                    raise RuntimeError,'Surface not supported yet'
            else:
                log.debug( '    %s.fireSingle Case 4: no surface nor neighbors nearby' % ('egfrd2' ) )
                # Case 4.
                pass

        # We can get here if:
        # 1. there were closeNeighbors. But it wasn't necessary to build a
        # multi (or just a pair was formed, todo), and there are bursted
        # singles to be restored.
        # 2. there was not a closeNeighbor nor a surface nearby.

        if distanceToSurface < closestShellDistance:
            closest = closestSurface
            closestShellDistance = distanceToSurface

        self.updateSingle( single, closest, closestShellDistance )

        bursted = uniq( bursted )
        burstedSingles = [ s for s in bursted if isinstance( s, Single ) ]
        # Probably everything would work just fine still if we didn't restore
        # here, since those burstedSingles are already in the scheduler with a
        # dt=0 because they have been given a shell with size minSize.
        self.restoreSingleShells( burstedSingles )
            
        log.info( 'single shell %g dt %g.' % ( single.size, single.dt ) )

        # Important. By returning 'dt', the eventScheduler reschedules this
        # single.
        return single.dt


    def createCylinder( self, single, pos, orientation, radius_a, r0, radius_b, halfLength, interactionType ):
        log.debug( '    %s.createCylinder: %s ' % ('egfrd2', single) )

        assert single.dt == 0.0
        assert single.getMobilitySize() == 0.0

        particle = single.particle
        rt = self.getReactionType1( particle.species )

        cylinder = Cylinder2( particle, rt, interactionType, pos, orientation, radius_a, r0, radius_b, halfLength )
        cylinder.initialize( self.t )

        return cylinder


    # Find largest possible cylinder around particle, such that it is not
    # interfering with other particles. Miedema's algorithm.
    def formCylinder( self, single, cylindricalSurface, interactionType ):
        log.debug( '    %s.formCylinder: %s' % ('egfrd2', single) )

        orientation = cylindricalSurface.orientation

        # Calculate origin (a.k.a. pos, or posPrime here) of cylinder.
        posPrime = cylindricalSurface.projection( single.pos )
        # Todo. Applyboundary.
        assert numpy.linalg.norm(single.pos - posPrime) >= cylindricalSurface.radius

        x0Prime = single.pos - posPrime
        z0 = numpy.dot( x0Prime, orientation )
        r0 = numpy.linalg.norm( x0Prime )

        dr = self.getMaxShellSize() #INF #CYLINDRICAL_SHELL_MAX_WIDTH - r0
        dz = self.getMaxShellSize() #INF #CYLINDRICAL_SHELL_MAX_HALF_LENGTH

        allNeighbors, distances = self.getNeighbors( single.pos, dr,
                                                       ignore=[single,] )
        closest = allNeighbors.pop()
        closestShellDistance = distances[-1]
        bursted = self.burstNonMultis( allNeighbors )
        for neighbor in bursted:
            rhoi = neighbor.shellList[0].size
            xiPrime = neighbor.shellList[0].origin - posPrime
            zi = numpy.dot( xiPrime, orientation )
            dzi = abs(zi) - rhoi
            temp = zi*numpy.array(orientation)
            temp2 = numpy.array(xiPrime) - numpy.array(temp)
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

        radius_a = cylindricalSurface.radius
        radius_b = r0 + dr
        halfLength = abs(dz)
        # Make cylinder with radius b and half-length dz. 
        cylinder = self.createCylinder( single, posPrime, orientation, radius_a, r0, radius_b, halfLength, interactionType )

        self.removeFromShellMatrix( single )

        self.shellMatrix.addCylinder( (cylinder, 0), cylinder.shellList[0] )
        cylinder.determineNextEvent( self.t )
        self.addCylinderEvent( cylinder )
        log.info( 'cylinder shell pos %s radius %g halfLength %g dt %g' % (cylinder.origin, cylinder.radius, cylinder.halfLength, cylinder.dt) )
        #Todo:
        #assert self.checkObj( cylinder )
        return cylinder
   

    def addCylinderEvent( self, cylinder ):
        print 'Cylinder!'
        log.debug( '    %s.addCylinderEvent: %s %g' % ('egfrd2', cylinder, cylinder.dt) )
        eventID = self.addEvent( self.t + cylinder.dt, 
                                 Delegate( self, EGFRDSimulator.fireCylinder ), 
                                 cylinder )
        cylinder.eventID = eventID


    def fireCylinder( self, cylinder ):
        log.debug( '    %s.fireCylinder: %s' % ('egfrd2', cylinder) )

        # Todo: is this ok here?
        self.shellMatrix.removeCylinder( cylinder )

        # Todo
        #assert self.checkObj( cylinder )

        log.info( 'fire: %s eventType %s' % (cylinder, cylinder.eventType ) )

        particle = cylinder.particle

        # Todo: needed?
        # self.applyBoundary( pos )

        # Three cases:
        #  Escaping (1) / interacting (0) through r.
        #  Escaping (2) through z.
        #  Single reaction (3)

        # 3. Single reaction.
        # Todo: copy-pasted from fireSingle, almost identical to firePair.
        # Move to new method.
        if cylinder.eventType == 3:

            log.info( 'single(cylinder) reaction %s' % str( cylinder ) )

            self.burstCylinder( cylinder )

            try:
                # Todo. Does this work:
                self.fireSingleReaction( cylinder )
            except NoSpace:
                log.info( 'single(reaction) reaction; placing product failed.' )
                # Todo. Do something like this.
                #self.addToShellMatrix( cylinder )
                #self.mainSimulator.rejectedMoves += 1
                #cylinder.reset()
                #return cylinder.dt
                raise Exception, 'Not implemented'

            cylinder.dt = -INF
            return cylinder.dt


        # 2. Escaping through z.
        elif cylinder.eventType == 2:

            log.info( 'cylinder %s escape through z' % (cylinder ) )
            #log.debug( 'r0 = %g, dt = %g, %s' %
            #               ( r0, cylinder.dt, cylinder.pgf.dump() ) )


            # calculate new r
            r = cylinder.drawR( cylinder.dt )

            # calculate displacement in theta
            dtheta = cylinder.drawTheta( cylinder.dt )

            # decide on displacement in z.
            if round(numpy.random.uniform()):
                dz = -cylinder.halfLength
            else:
                dz = cylinder.halfLength


        # 1 or 0. Escaping (1) / interacting (0) through r.
        elif cylinder.eventType == 0:

            # calculate displacement in z
            dz = cylinder.drawZ( cylinder.dt )

            # calculate new r
            cylinder.eventType = cylinder.drawEventType( cylinder.dt )
            if cylinder.eventType == 0:
                # Interaction
                log.info( 'cylinder %s interaction with surface' % ( cylinder ) )
                try:
                    # Todo. Does this work:
                    self.fireSingleReaction( cylinder, cylinder.interactionType )
                except NoSpace:
                    log.info( 'single(reaction) reaction; placing product failed.' )
                    raise Exception, 'Not implemented'

                cylinder.dt = -INF
                return cylinder.dt
            else:
                # Escape
                log.info( 'cylinder %s escape through r' % ( cylinder ) )
                r = cylinder.radius_b

                # calculate displacement in theta
                dtheta = cylinder.drawTheta( cylinder.dt )

        else:
            raise SystemError, 'Bug: invalid eventType.'


        # calculate total displacement
        newpos = cylinder.newPosition( r, dtheta, dz )
        self.applyBoundary( newpos )

        assert self.checkOverlap( newpos, particle.species.radius,
                                  ignore = [ particle ] )

        newsingle = self.createSingle( particle )
        self.addToShellMatrix( newsingle )
        self.addSingleEvent( newsingle )

        assert self.checkObj( newsingle )

        cylinder.dt = -INF
        return cylinder.dt

    """
