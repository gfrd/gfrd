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
    # Todo. Should also depend on reaction rates. dt * k << 1.
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
    log.debug( '\t\tDebug. bd dt = %.3g' % dt )

    return dt


class BDSimulatorCoreBase( object ):
    """BDSimulatorCoreBase.

    Child classes of BDSimulatorCoreBase are:
    * BDSimulatorCore
        - Instantiated by BDSimulator as BDSimulator.core.
        - Uses the particleMatrix from the simulator that instantiates it, 
          that would be a BDSimulator, for particle detection. BDsimulator 
          derives from ParticleSimulatorBase, which has declared this 
          particleMatrix.
    * MultiBDCore
        - Instantiated by Multi as Multi.core.
        - Uses a self declared particleMatrix for particle detection, which 
          contains only the particles that are in the Multi, not other 
          particles that might still be in the EGFRDSimulator.

    BDSimulatorCoreBase uses a particle*List* to loop over all particles in 
    each step, and propagate them.

    BDSimulatorCoreBase borrows the following from the main simulator:
    - speciesList
    - reactionTypes list, both 1 and 2.
    
    """

    def __init__( self, main ):
        # Reference to the main, egfrd, simulator.
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


    def getP_acct( self, rt, D, sigma ):
        # Todo. This does not obey detailed balance yet in 1D and 2D and with 
        # surfaces.
        try:
            return self.P_acct[rt]

        except KeyError:
            I = _gfrd.I_bd( sigma, self.dt, D )
            p = rt.k * self.dt / ( I * 4.0 * numpy.pi )
            if not 0.0 <= p < 1.0:
                log.error( 'Invalid acceptance ratio (%s) for '
                           'reaction %s. rt.k = %.3g, I = %.3g' %
                           ( p, rt, rt.k, I ) )
                p = 1.0
            self.P_acct[rt] = p
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
            # Todo. What if particle has been removed after a reaction?
            # Check if particle still exists in original particleList or
            # keep exclude list.
            particle = self.particlesToStep.pop() # take the last one
            self.propagateParticle( particle )


    def propagateParticle( self, particle ):
        species = particle.species

        # 1. Try single reactions first.
        rt1 = self.attemptSingleReactions( species )
        if rt1:
            # Todo. Shouldn't we displace the particle first? 
            try:
                self.fireReaction1( particle, rt1 )
            except NoSpace:
                log.info( '\t\tfireReaction1 rejected.' )
            return

        D = species.D
        if D == 0.0:
            return

        displacement = particle.surface.drawBDdisplacement( self.dt, D )
        log.debug( '\t\tDebug. Multi. %s, displacement = [%.3g, %.3g, %.3g].' %
                   ( particle, displacement[0], displacement[1], 
                     displacement[2] ) )

        newpos = particle.pos + displacement
        newpos %= self.main.worldSize   #self.applyBoundary( newpos )
        
        neighbors = self.getParticlesWithinRadiusNoSort( newpos, species.radius,
                                                         ignore=[ particle, ] )

        # 2. Try reaction of 2 particles, even if newpos also overlaps with 
        # surface. 
        if neighbors:
            if len( neighbors ) >= 2:
                log.info( '\t\tcollision two or more particles; move rejected' )
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
                    log.info( '\t\tfire reaction2' )
                    try:
                        self.fireReaction2( particle, closest, rt )
                    except NoSpace:
                        log.info( '\t\tfireReaction2 move rejected' )
                    return

            else:
                log.info( '\tcollision move rejected' )

            # Neighbor is reflecting us. Don't move.
            return

        # 3. Try binding with surface.
        distanceToSurface, surface = \
            self.main.getClosestSurfaceWithinRadius( newpos, species.radius,
                                                     ignore=[ particle, ] )
        if surface:
            rt = self.main.getInteractionType( species, surface )

            if rt.k != 0.0:
                radius12 = species.radius + surface.Lz
                p = self.getP_acct( rt, D, radius12 )

                rnd = numpy.random.uniform()

                if p > rnd:
                    log.info( '\tfire interaction' )
                    try:
                        self.fireReaction1( particle, rt )
                    except NoSpace:
                        log.info( '\tfireReaction1 move rejected' )
                    return

            else:
                log.info( '\tinteraction move rejected' )

            # Surface is reflecting us. Don't move.
            return



        # No reaction or interaction. 
        try:
            self.clearVolume( newpos, particle.radius, ignore=[ particle, ] )
            self.moveParticle( particle, newpos )
        except NoSpace:
            log.info( '\tpropagation move rejected.' )


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

            self.lastReaction = Reaction( rt, [ particle ], [] )
            
        elif len( rt.products ) == 1:
            
            productSpecies = rt.products[0]
            radius = productSpecies.radius
    
            if( self.main.isSurfaceUnbindingReaction( rt, currentSurface ) ):
                newpos = currentSurface.randomUnbindingSite( oldpos, radius )

            elif( self.main.isSurfaceBindingReaction( rt, currentSurface ) ):
                # Select position on surface with z = 0.
                newpos, _ = productSpecies.surface.projectedPoint( oldpos )
            else:
                newpos = oldpos

            if not self.checkOverlap( newpos, radius,
                                      ignore=[ particle, ] ):
                log.info( '\tno space for product particle.' )
                raise NoSpace()

            self.clearVolume( newpos, radius, ignore=[ particle, ] )
                
            self.removeParticle( particle )
            newparticle = self.createParticle( productSpecies, newpos )

            self.lastReaction = Reaction( rt, [ particle ], [ newparticle ] )

            
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

                if( self.main.isDirectSurfaceUnbindingReaction( rt ) ):
                    newpos1 = oldpos
                    # Particle2 ends up in world (defaultSurface).
                    newpos2 = currentSurface.randomUnbindingSite( oldpos, 
                                                                  pairDistance )
                else:
                    # (1.0 + 1e-10) # safety
                    vector = surface.randomVector( pairDistance )
                
                    # place particles according to the ratio D1:D2
                    # this way, species with D == 0 doesn't move.
                    # FIXME: what if D1 == D2 == 0?
                    newpos1 = oldpos + vector * ( D1 / D12 )
                    newpos2 = oldpos - vector * ( D2 / D12 )
                    #FIXME: check surfaces here

            
                self.applyBoundary( newpos1 )
                self.applyBoundary( newpos2 )

                # accept the new positions if there is enough space.
                if self.checkOverlap( newpos1, radius1, 
                                      ignore=[ particle, ] ) and \
                   self.checkOverlap( newpos2, radius2, ignore=[ particle, ] ):
                    break
            else:
                log.info( '\tno space for product particles.' )
                raise NoSpace()

            self.clearVolume( newpos1, radius1, ignore=[ particle, ] )
            self.clearVolume( newpos2, radius2, ignore=[ particle, ] )

            # move accepted
            self.removeParticle( particle )

            newparticle1 = self.createParticle( productSpecies1, newpos1 )
            newparticle2 = self.createParticle( productSpecies2, newpos2 )

            self.lastReaction = Reaction( rt, [ particle ], 
                                          [ newparticle1, newparticle2 ] )

        else:
            raise RuntimeError( 'num products >= 3 not supported.' )

        self.reactionEvents += 1


    def fireReaction2( self, particle1, particle2, rt ):
        #if( self.main.isDirectSurfaceBindingReaction( rt ) ):
            # Todo. Direct binding.

        assert particle1.surface == particle2.surface

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
            newparticle = self.createParticle( productSpecies, newPos )

            try:
                self.particlesToStep.remove( particle2 )
            except ValueError:  
                pass     # particle2 already stepped, which is fine.

            self.reactionEvents += 1

            self.lastReaction = Reaction( rt, [ particle1, particle2 ], 
                                          [ newparticle ] )

            return
        
        else:
            raise NotImplementedError( 'num products >= 2 not supported.' )


    def check( self ):
        # particles don't overlap
        for particle in self.particleList:
            assert self.checkOverlap( particle.pos, particle.radius,
                                      ignore=[ particle, ] )



##########################################################################
class BDSimulatorCore( BDSimulatorCoreBase ):
    """Used by BDSimulator only.  

    """

    def __init__( self, main ):
        # main can be a BDSimulator object. 
        BDSimulatorCoreBase.__init__( self, main )

        self.checkOverlap = self.main.checkOverlap

        self.moveParticle = self.main.moveParticle

        #self.getNeighborParticles = main.getNeighborParticles
        self.getParticlesWithinRadius = main.getParticlesWithinRadius
        self.getParticlesWithinRadiusNoSort = \
            main.getParticlesWithinRadiusNoSort
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


    def createParticle( self, species, pos ):
        particle = self.main.createParticle( species, pos )
        self.addToParticleList( particle )


    def clearVolume( self, pos, radius, ignore=[] ):
        """This method is a customization point for implementing BD in 
        protective domains.

        """
        pass


class BDSimulator( ParticleSimulatorBase ):
    """Can be used for a Brownian Dynamics simulation.

    """

    def __init__( self, worldSize ):
        ParticleSimulatorBase.__init__( self, worldSize )

        self.core = BDSimulatorCore( self )
        self.isDirty = True


    def gett( self ):
        return self.core.t
    def sett( self, t ):
        self.core.t = t
    t = property( gett, sett )


    def getDt( self ):
        return self.core.dt
    dt = property( getDt )


    def getStepCounter( self ):
        return self.core.stepCounter
    stepCounter = property( getStepCounter )


    def initialize( self ):
        self.setAllRepulsive()

        self.core.initialize()

        self.isDirty = False


    def getNextTime( self ):
        return self.core.t + self.core.dt


    def reset( self ):
        # DUMMY
        self.core.t = 0


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

        log.info( '\t%d: t = %.3g dt = %.3g, reactions = %d, rejectedMoves = %d'
                  % ( self.stepCounter, self.t, self.dt, 
                      self.reactionEvents, self.rejectedMoves ) )


    def check( self ):
        pass
