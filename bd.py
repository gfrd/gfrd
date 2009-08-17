#!/usr/env python

import weakref

import math

import numpy


from utils import *
from surface import *

from gfrdbase import *
import _gfrd

DEFAULT_DT_FACTOR = 1e-5

def calculateBDDt( speciesList, factor ):

    D_list = []
    radius_list = []
    for species in speciesList:
        if species.pool.size != 0:
            D_list.append( species.D )
            radius_list.append( species.radius )
    D_max = max( D_list ) * 2  # max relative diffusion speed
    radius_min = min( radius_list )
    sigma_min = radius_min * 2

    dt = factor * sigma_min ** 2 / D_max  
    log.debug( 'bd dt = %g' % dt )

    return dt

'''
Used by 
- BDSimulatorCore (for BDSimulator).
- MultiBDCore (for Multi.core)
'''
class BDSimulatorCoreBase( object ):
    
    '''
    BDSimulatorCoreBase borrows the following from the main simulator:
    - speciesList
    - reactionTypes list (both 1 and 2)
    
    '''
    def __init__( self, main ):

        # Reference to the main (egfrd) simulator.
        self.main = weakref.proxy( main )

        self.particleList = []

        self.t = 0.0
        self.dt = 0.0

        self.dtFactor = DEFAULT_DT_FACTOR

        self.stepCounter = 0

        self.reactionEvents = 0

        self.speciesList = self.main.speciesList
        self.getReactionType1 = self.main.getReactionType1
        self.getReactionType2 = self.main.getReactionType2
        self.applyBoundary = self.main.applyBoundary

        self.lastReaction = None

        self.P_acct = {}


    def initialize( self ):
        self.determineDt()

    def clearParticleList( self ):
        self.particleList = []

    def addToParticleList( self, particle ):
        self.particleList.append( particle )

    def removeFromParticleList( self, particle ):
        self.particleList.remove( particle )


    def getNextTime( self ):
        return self.t + self.dt

    def stop( self, t ):
        # dummy
        self.t = t

    def determineDt( self ):

        self.dt = calculateBDDt( self.speciesList.values(), self.dtFactor )


    # Todo. Does this obey detailed balance?
    def getP_acct( self, rt, D, sigma ):

        try:
            return self.P_acct[ rt ]

        except KeyError:
            I = _gfrd.I_bd( sigma, self.dt, D )
            p = rt.k * self.dt / ( I * 4.0 * numpy.pi )
            if not 0.0 <= p < 1.0:
                raise RuntimeError,\
                    'Invalid acceptance ratio (%s) for reaction %s.' \
                    % ( p, rt )
            self.P_acct[ rt ] = p
            return p


    def step( self ):

        self.stepCounter += 1
        self.lastReaction = None

        self.propagate()

        self.t += self.dt


    def propagate( self ):
        
        # Copy particlelist, because it may change
        self.particlesToStep = self.particleList[:]

        random.shuffle( self.particlesToStep )
        while self.particlesToStep:
            particle = self.particlesToStep.pop() # take the last one
            self.propagateParticle( particle )


    def propagateParticle( self, particle ):

        species = particle.species

        rt1 = self.attemptSingleReactions( species )
        if rt1:
            # Todo. Shouldn't we propagate first???
            try:
                self.fireReaction1( particle, rt1 )
            except NoSpace:
                log.info( 'fireReaction1 rejected.' )
            return

        D = species.D
        if D == 0.0:
            return

        displacement = particle.surface.drawBDdisplacement( self.dt, D )

        newpos = particle.pos + displacement
        newpos %= self.main.worldSize   #self.applyBoundary( newpos )
        
        neighbors = self.getParticlesWithinRadiusNoSort( newpos, species.radius,
                                                         ignore=[particle] )
        '''
        Try reaction of 2 particles first, even if newpos also overlaps with surface. 
        '''
        if neighbors:

            if len( neighbors ) >= 2:
                log.info( 'collision two or more particles; move rejected' )
                return

            closest = neighbors[0]
            species2 = closest.species

            rt = self.main.reactionTypeMap2.get( ( species, species2 ) )

            if rt.k != 0.0:
                radius12 = species.radius + species2.radius
                D12 = D + species2.D

                p = self.getP_acct( rt, D12, radius12 )

                rnd = numpy.random.uniform()

                if p > rnd:
                    log.info( 'fire reaction2' )
                    try:
                        self.fireReaction2( particle, closest, rt )
                    except NoSpace:
                        log.info( 'fireReaction2 move rejected' )
                    return

            else:
                log.info( 'collision move rejected' )

            # Neighbor is reflecting us. Don't move.
            return

        surface, distanceToSurface = self.main.getClosestSurfaceWithinRadius( newpos, species.radius, ignore=[particle] )
        if surface:
            rt = self.main.getInteractionType( species, surface )

            if rt.k != 0.0:
                radius12 = species.radius + surface.Lz
                p = self.getP_acct( rt, D, radius12 )

                rnd = numpy.random.uniform()

                if p > rnd:
                    log.info( 'fire interaction' )
                    try:
                        self.fireReaction1( particle, rt )
                    except NoSpace:
                        log.info( 'fireReaction1 move rejected' )
                    return

            else:
                log.info( 'interaction move rejected' )

            # Surface is reflecting us. Don't move.
            return




        # No reaction or interaction. 
        try:
            self.clearVolume( newpos, particle.radius, ignore=[particle] )
            self.moveParticle( particle, newpos )
        except NoSpace:
            log.info( 'propagation move rejected.' )



    def attemptSingleReactions( self, species ):

        reactionTypes = self.getReactionType1( species )
        if not reactionTypes:
            return None  # no reaction

        rnd = numpy.random.uniform() / self.dt

        # handle the most common case efficiently.
        if len( reactionTypes ) == 1:  
            if reactionTypes[0].k >= rnd:
                return reactionTypes[0]
            else:
                return None

        # if there are more than one possible reaction types..
        k_array = numpy.add.accumulate( [ rt.k for rt in reactionTypes ] )

        if k_array[-1] < rnd:
            return None

        i = numpy.searchsorted( k_array, rnd )

        return reactionTypes[i]


    def fireReaction1( self, particle, rt ):
        currentSurface = particle.surface
        
        oldpos = particle.pos.copy()

        if len( rt.products ) == 0:
            
            self.removeParticle( particle )

            self.lastReaction = Reaction( rt, [particle], [] )
            
        elif len( rt.products ) == 1:
            
            productSpecies = rt.products[0]
            radius = productSpecies.radius
    
	    if isinstance( rt, SurfaceUnbindingReactionType ):
		newSurface = self.main.defaultSurface
		newpos = currentSurface.randomUnbindingSite( oldpos, radius )
            elif isinstance( rt, SurfaceBindingInteractionType ):
                # Get interaction surface from reactants list.
                newSurface = rt.reactants[1]
                # Todo. Does this obey detailed balance?
                # Select position on surface with z=0.
                newpos, _ = newSurface.calculateProjection( oldpos )
            else:
                newSurface = currentSurface
                newpos = oldpos

            if not self.checkOverlap( newpos, radius,
                                      ignore = [ particle, ] ):
                log.info( 'no space for product particle.' )
                raise NoSpace()

            self.clearVolume( newpos, radius, ignore = [ particle ] )
                
            self.removeParticle( particle )
            newparticle = self.createParticle( productSpecies, newpos, newSurface )

            self.lastReaction = Reaction( rt, [particle], [newparticle] )

            
        elif len( rt.products ) == 2:
            
            productSpecies1 = rt.products[0]
            productSpecies2 = rt.products[1]
            
            D1 = productSpecies1.D
            D2 = productSpecies2.D
            D12 = D1 + D2
            
            radius1 = productSpecies1.radius
            radius2 = productSpecies2.radius
            radius12 = radius1 + radius2

            for i in range( 100 ):
                rnd = numpy.random.uniform()
                pairDistance = drawR_gbd( rnd, radius12, self.dt, D12 )

                if isinstance( rt, SurfaceDirectUnbindingReactionType ):
                    # Direct unbinding.
                    # Particle2 ends up in world (defaultSurface).
                    newSurface1 = currentSurface
                    newSurface2 = self.main.defaultSurface
                    newpos1 = oldpos
                    # Todo. Does this obey detailed balance?
                    newpos2 = currentSurface.randomUnbindingSite( oldpos, pairDistance )
                else:
                    vector = surface.randomVector( pairDistance ) # (1.0 + 1e-10) # safety
                
                    # place particles according to the ratio D1:D2
                    # this way, species with D=0 doesn't move.
                    # FIXME: what if D1 == D2 == 0?
                    newpos1 = oldpos + vector * ( D1 / D12 )
                    newpos2 = oldpos - vector * ( D2 / D12 )
                    newSurface1 = currentSurface
                    newSurface2 = currentSurface
                    #FIXME: check surfaces here

            
                self.applyBoundary( newpos1 )
                self.applyBoundary( newpos2 )

                # accept the new positions if there is enough space.
                if self.checkOverlap( newpos1, radius1,
                                      ignore = [ particle, ]) and \
                                      self.checkOverlap( newpos2, radius2,
                                                         ignore = 
                                                         [ particle, ]):
                    break
            else:
                log.info( 'no space for product particles.' )
                raise NoSpace()

            self.clearVolume( newpos1, radius1, ignore = [ particle ] )
            self.clearVolume( newpos2, radius2, ignore = [ particle ] )

            # move accepted
            self.removeParticle( particle )

            newparticle1 = self.createParticle( productSpecies1, newpos1, newSurface1 )
            newparticle2 = self.createParticle( productSpecies2, newpos2, newSurface2 )

            self.lastReaction = Reaction( rt, [particle], 
                                          [newparticle1, newparticle2] )

        else:
            raise RuntimeError, 'num products >= 3 not supported.'

        self.reactionEvents += 1




    def fireReaction2( self, particle1, particle2, rt ):
        # Todo. Direct binding.
        assert particle1.surface == particle2.surface
        newSurface = particle1.surface

        pos1 = particle1.pos.copy()
        pos2 = particle2.pos.copy()

        if len( rt.products ) == 1:
                
            productSpecies = rt.products[0]

            D1 = particle1.species.D
            D2 = particle2.species.D

            pos2t = cyclicTranspose( pos2, pos1, self.main.worldSize )
            newPos = ( D2 * pos1 + D1 * pos2t ) / ( D1 + D2 )
            self.applyBoundary( newPos )

            if not self.checkOverlap( newPos, productSpecies.radius,
                                      ignore=[ particle1, particle2 ] ):
                raise NoSpace()
            self.clearVolume( newPos, productSpecies.radius,
                              ignore=[ particle1, particle2 ] )

            self.removeParticle( particle1 )
            self.removeParticle( particle2 )
            newparticle = self.createParticle( productSpecies, newPos, newSurface )

            try:
                self.particlesToStep.remove( particle2 )
            except ValueError:  
                pass     # particle2 already stepped, which is fine.

            self.reactionEvents += 1

            self.lastReaction = Reaction( rt, [particle1, particle2], 
                                          [newparticle] )

            return
        
        else:
            raise NotImplementedError,\
                'num products >= 2 not supported.'


    def check( self ):

        # particles don't overlap

        for particle in self.particleList:
            assert self.checkOverlap( particle.pos, particle.radius,
                                      ignore=[particle,] )



##########################################################################
'''
Used by BDSimulator only.  
'''
class BDSimulatorCore( BDSimulatorCoreBase ):

    def __init__( self, main ):

        # main can be a BDSimulator object. 
        BDSimulatorCoreBase.__init__( self, main )

        self.checkOverlap = self.main.checkOverlap

        self.moveParticle = self.main.moveParticle

        #self.getNeighborParticles = main.getNeighborParticles
        self.getParticlesWithinRadius = main.getParticlesWithinRadius
        self.getParticlesWithinRadiusNoSort = main.getParticlesWithinRadiusNoSort
        #self.getClosestParticle = main.getClosestParticle

        
    def initialize( self ):
        BDSimulatorCoreBase.initialize( self )

        # BDSimulator uses a particle list instead of the particleMatrix.
        self.updateParticleList()

    def updateParticleList( self ):
        self.clearParticleList()

        particles = self.main.particleMatrix.getAll( )
        # This way particle.surface is also known.
        for particle in particles:
            self.addToParticleList( particle )

    def addParticle( self, particle ):
        self.main.addParticle( particle )
        self.addToParticleList( particle )

    def removeParticle( self, particle ):
        self.main.removeParticle( particle )
        self.removeFromParticleList( particle )

    def createParticle( self, species, pos, surface ):
        particle = self.main.createParticle( species, pos, surface )
        self.addToParticleList( particle )


    # This method is a customization point for implementing
    # BD in protective domains.

    def clearVolume( self, pos, radius, ignore=[] ):
        
        pass


'''
Can be used for a Brownian Dynamics simulation.
'''
class BDSimulator( ParticleSimulatorBase ):
    
    def __init__( self ):

        ParticleSimulatorBase.__init__( self )

        self.core = BDSimulatorCore( self )
        self.isDirty = True

    def gett( self ):
        return self.core.t

    def sett( self, t ):
        self.core.t = t

    def getDt( self ):
        return self.core.dt

    def getStepCounter( self ):
        return self.core.stepCounter

    t = property( gett, sett )
    dt = property( getDt )
    stepCounter = property( getStepCounter )

    def initialize( self ):

        self.setAllRepulsive()

        self.core.initialize()

        self.isDirty = False


    def getNextTime( self ):
        return self.core.t + self.core.dt

    def reset( self ):
        # DUMMY
        self.core.t=0

    def stop( self, t ):
        # dummy
        self.core.stop( t )

    def step( self ):

        self.reactionType = None

        if self.isDirty:
            self.initialize()

        #if self.stepCounter % 10000 == 0:
        #    self.check()

        self.core.step()

        log.info( '%d: t=%g dt=%g, reactions=%d, rejectedMoves=%d' %
                  ( self.stepCounter, self.t, self.dt, self.reactionEvents,
                    self.rejectedMoves ) )

    def check( self ):
        pass
