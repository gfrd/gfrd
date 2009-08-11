#!/usr/env python


import math
import random

import numpy
import scipy


from utils import *
from surface import *
# _gfrd is a library that is written in C++ and pythonified using Boost 
# library.
from _gfrd import *

from cObjectMatrix import *

from log import *

def p_free( r, t, D ):
    Dt4 = D * t * 4.0
    Pi4Dt = numpy.pi * Dt4
    rsq = r * r
    
    p = math.exp( - rsq / Dt4 ) / math.sqrt( Pi4Dt * Pi4Dt * Pi4Dt )

    jacobian = 4.0 * numpy.pi * rsq

    return p * jacobian
    

class NoSpace( Exception ):
    pass


class Species( object ):
    
    def __init__( self, id, D, radius ):
        self.id = id
        self.D = D
        self.radius = radius
        self.pool = ParticlePool() # Stores positions only.


    def newParticle( self, position ):

        serial = self.pool.newParticle( position )
        return serial


    def removeParticleByIndex( self, index ):
        self.pool.removeByIndex( index )

    def removeParticleBySerial( self, serial ):
        self.pool.removeBySerial( serial )


class ReactionType( object ):

    def __init__( self, reactants=[], products=[], k=0.0 ):
        self.reactants = reactants
        self.products = products
        self.k = k
        self.check()


    def check( self ):
        if len( self.reactants ) == 1 and len( self.products ) <= 2:
            totalProductRadii = 0.0
            for product in self.products:
                totalProductRadii += product.radius
            if totalProductRadii > self.reactants[0].radius * 2:
                raise ValueError, \
                    'total product radii must be smaller than ' \
                    + 'reactant radius * 2'
        if self.k < 0:
            raise ValueError, 'k < 0'

    def order( self ):
        return len( self.reactants )

    def __str__( self ):
        s = ''
        for i in self.reactants:
            if isinstance( i, Surface ):
                s += str( i )
            else:
                s += i.id
            s += ' '

        s += '-> '

        for i in self.products:
            s += i.id
            s += ' '

        return s[:-1]


'''
Reaction Types.

A -> B
'''
class UnimolecularReactionType( ReactionType ):

    def __init__( self, s1, p1, k ):
        ReactionType.__init__( self, [ s1, ], [ p1, ], k )


'''
A -> None
'''
class DecayReactionType( ReactionType ):

    def __init__( self, s1, k ):
        ReactionType.__init__( self, [ s1, ], [], k )


'''
A + B -> C
'''
class BindingReactionType( ReactionType ):

    def __init__( self, s1, s2, p1, k ):
        ReactionType.__init__( self, [ s1, s2 ], [ p1, ], k )
        # Todo. These are not used.
        D = s1.D + s2.D
        sigma = s1.radius + s2.radius


'''
C -> A + B
'''
class UnbindingReactionType( ReactionType ):

    def __init__( self, s1, p1, p2, k ):
        ReactionType.__init__( self, [ s1, ], [ p1, p2 ], k )
 

'''
A + B is repulsive.
'''
class RepulsionReactionType( ReactionType ):

    def __init__( self, s1, s2 ):
        ReactionType.__init__( self, [ s1, s2 ], [], 0.0 )
        # Todo. These are not used.
        D = s1.D + s2.D
        sigma = s1.radius + s2.radius


'''
Surface binding.

Surface binding is not a Poisson process, so interactions should not be added 
to the reaction list using addReactionType but added to the interaction list using
addInteractionType.

A + Surface -> B_on_Surface
'''
class SurfaceBindingInteractionType( ReactionType ):
    def __init__( self, reactantSpecies, surface, productSpecies,  k ):
        ReactionType.__init__( self, [ reactantSpecies, surface ], [ productSpecies, ], k )


'''
A + B_on_Surface + Surface -> C_on_Surface
A + Surface should be repulsive.
'''
class SurfaceDirectBindingInteractionType( ReactionType ):
    def __init__( self, reactantSpecies1, reactantSpecies2, surface, productSpecies,  k ):
        ReactionType.__init__( self, [ reactantSpecies1, reactantSpecies2, surface ], [ productSpecies, ], k )


'''
A + Surface is repulsive.
'''
class SurfaceRepulsionInteractionType( ReactionType ):
    def __init__( self, species, surface ):
        ReactionType.__init__( self, [ species, surface ], [], 0.0 )


'''
Surface unbinding.

Surface unbinding is a Poisson process, so these reactions can be added to the 
reaction list by using addReactionType.

A_on_Surface -> B
'''
class SurfaceUnbindingReactionType( ReactionType ):
    def __init__( self, reactantSpecies, productSpecies, k ):
        # After unbinding from a surface the particle always ends up on the 
        # defaultSurface (world) for now.
        ReactionType.__init__( self, [ reactantSpecies ], [ productSpecies, ], k )


'''
A_on_Surface -> B_on_Surface + C
'''
class SurfaceDirectUnbindingReactionType( ReactionType ):
    def __init__( self, reactantSpecies, productSpecies1, productSpecies2, k ):
        # After unbinding from a surface, particle2 always ends up on the 
        # defaultSurface (world), and particle1 stays on the surface it was 
        # on.
        ReactionType.__init__( self, [ reactantSpecies ], [ productSpecies1, productSpecies2 ], k )


##############################################################################

class Reaction:
    def __init__( self, type, reactants, products ):
        self.type = type
        self.reactants = reactants
        self.products = products

    def __str__( self ):
        return 'Reaction( ' + str( self.type ) + ', ' + str( self.reactants )\
            + ', ' + str( self.products ) + ' )'


class Particle( object ):

    ### Not only constructor, but can also be used to make a copy of a
    ### particle by index or serial.
    def __init__( self, species, serial=None, index=None, surface=None ):

        self.species = species

        if not serial is None:
            self.serial = serial
        elif not index is None:
            self.serial = species.pool.getSerialByIndex( index )
        else:
            raise ValueError, 'give either serial or index.'

        self.hash = hash( self.species ) ^ self.serial

        self.surface = surface

    def __str__( self ):

        return "( '" + self.species.id + "', " + str( self.serial ) + ' )'

    def __repr__( self ):

        return self.__str__()


    def __eq__( self, other ):

        return self.species == other.species and self.serial == other.serial

    def __ne__( self, other ):

        return self.species != other.species or self.serial != other.serial

    def __cmp__( self, other ):

        if self.species == other.species:
            return self.serial - other.serial
        elif self.species < other.species:
            return -1
        else:
            return 1


    def __hash__( self ):

        return self.hash


    def getPos( self ):
        pool = self.species.pool
        return pool.positions[ pool.indexMap[ self.serial ] ]

    def setPos( self, newpos ):
        pool = self.species.pool
        pool.positions[ pool.indexMap[ self.serial ] ] = newpos

    pos = property( getPos, setPos )

    def getRadius( self ):
        return self.species.radius

    radius = property( getRadius )


    def getIndex( self ):
        return self.species.pool.indexMap[ self.serial ]


class DummyParticle( object ):
    def __init__( self ):
        self.species = None
        self.serial = -1



# Stores positions only.
class ParticlePool( object ):

    def __init__( self ):

        self.indexMap = {}
        self.serialCounter = 0

        self.serials = numpy.array( [], numpy.integer )

        self.positions = numpy.array( [], numpy.floating )
        self.positions.shape = ( 0, 3 )

        self.size = 0


    def newParticle( self, position ):

        newindex = self.size
        newserial = self.serialCounter
        
        self.size += 1
        self.serialCounter += 1

        self.__resizeArrays( self.size )

        self.indexMap[ newserial ] = newindex
        self.serials[ newindex ] = newserial
        self.positions[ newindex ] = position

        return newserial
    

    def __resizeArrays( self, newsize ):

        self.serials.resize( newsize )
        self.positions.resize( ( newsize, 3 ) )


    def removeBySerial( self, serial ):

        index = self.getIndex( serial )
        self.removeByIndex( index )


    def removeByIndex( self, index ):

        self.size -= 1

        if self.size == 0:
            self.indexMap = {}
            self.__resizeArrays( self.size )
            return

        serialOfRemovedItem = self.serials[ index ]
        serialOfMovedItem = self.serials[ self.size ]

        # swap items at index and the end, and discard the item
        # to be removed, which is now at the end, by shrinking
        # the arrays by one.
        self.serials[ index ] = self.serials[ self.size ]
        self.positions[ index ] = self.positions[ self.size ]
        self.__resizeArrays( self.size )

        # book keeping
        del self.indexMap[ serialOfRemovedItem ]
        self.indexMap[ serialOfMovedItem ] = index



    def getSerialByIndex( self, index ):
        
        return self.serials[index]

    def getIndex( self, serial ):

        return self.indexMap[ serial ]



class ParticleSimulatorBase( object ):
    
    def __init__( self ):

        self.speciesList = {}
        ### Per species key a list of reactiontypes.
        self.reactionTypeMap1 = {}
        self.reactionTypeMap2 = {}
        self.interactionTypeMap = {}

        self.surfaceList = []

        #self.dt = 1e-7
        #self.t = 0.0

        self.H = 3.0
        
        self.dtLimit = 1e-3
        self.dtMax = self.dtLimit

        # counters
        self.rejectedMoves = 0
        self.reactionEvents = 0

        self.particleMatrix = SphereMatrix() # Stores complete particles (unlike particlePool).

        self.setWorldSize( INF )

        self.lastReaction = None

        # Needed when creating a shell around a particle whose surface is not 
        # specified. Backward compatibility.
        self.defaultSurface = World( )


    def initialize( self ):
        pass


    '''
    Returns sorted list of pairs:
    - distance to surface
    - surface itself

    We can not use objectmatrix, it would miss a surface if the origin of the 
    surface would not be in the same or neighboring cells as pos.
    '''
    def getClosestSurface( self, pos, ignore ):
        surfaces = [ None ]
        distances = [ INF ]
        ignoreSurfaces = []
        for obj in ignore:
            if isinstance( obj.surface, Surface ):
                # Ignore surface that particle is currently on.
                ignoreSurfaces.append( obj.surface )
	for surface in self.surfaceList:
            if surface not in ignoreSurfaces:
		posTransposed = cyclicTranspose( pos, surface.origin, self.worldSize )
	        distanceToSurface = surface.signedDistanceTo( posTransposed )
                if distanceToSurface < 0.0:
                    self.errors += 1
                distances.append( distanceToSurface )
		surfaces.append( surface )
        return min( zip( distances, surfaces ))


    def getClosestSurfaceWithinRadius( self, pos, radius, ignore ):
        distanceToSurface, closestSurface = self.getClosestSurface( pos, ignore ) 
        if distanceToSurface < radius:
            return  closestSurface, distanceToSurface
        else:
            return None, INF


    def reconstructParticleMatrix( self ):

        self.particleMatrix.clear()
        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                particle = Particle( species, index=i )
                self.particleMatrix.add( particle,
                                         particle.pos, species.radius )



    def setWorldSize( self, size ):

        if isinstance( size, list ) or isinstance( size, tuple ):
            size = numpy.array( size )

        self.worldSize = size

        self.particleMatrix.setWorldSize( size )

        # Note: the underscored functions here do not map directly to the 
        # underscored functions in distance.hpp.
        # distance_Simple.. -> utils.py -> pyGFRD.cpp -> distance.hpp
        # distance_Cyclic.. -> utils.py -> pyGFRD.cpp -> distance.hpp
        if isinstance( size, float ) and size == INF:
            log.debug( 'distance_Simple.. is used' )
            self._distance = distance_Simple
            #self._distanceArray = distanceArray_Simple
            self._distanceSq = distanceSq_Simple
            self._distanceSqArray = distanceSqArray_Simple
        else:
            log.debug( 'distance_Cyclic.. is used' )
            self._distance = distance_Cyclic
            #self._distanceArray = distanceSqArray_Cyclic
            self._distanceSq = distanceSq_Cyclic
            self._distanceSqArray = distanceSqArray_Cyclic

    def getWorldSize( self ):
        return self.worldSize

    def setMatrixSize( self, size ):
        self.particleMatrix.setMatrixSize( size )

    # Use this method to account for periodic boundary conditions.
    def applyBoundary( self, pos ):
        pos %= self.worldSize

    def getReactionType1( self, species ):
        return self.reactionTypeMap1.get( species, None )

    def getReactionType2( self, species1, species2 ):
        return self.reactionTypeMap2.get( ( species1, species2 ), None )

    def getInteractionType( self, species, surface ):
        return self.interactionTypeMap.get( (species, surface) )
    
    def getSpeciesByIndex( self, i ):
        return self.speciesList.values()[i]


    def distanceSq( self, position1, position2 ):
        return self._distanceSq( position1, position2, self.worldSize )

    def distance( self, position1, position2 ):
        return self._distance( position1, position2, self.worldSize )
        
    def distanceSqArray( self, position1, positions ):
        return self._distanceSqArray( position1, positions, self.worldSize )

    def distanceArray( self, position1, positions ):
        return numpy.sqrt( self.distanceSqArray( position1,\
                                                 positions ) )

    def addSurface( self, surface ):
        self.surfaceList.append( surface )


    def addSpecies( self, species ):
        self.speciesList[ species.id ] = species

    def addReactionType( self, rt ):

        numReactants = len( rt.reactants )

        if numReactants == 1:
            species1 = rt.reactants[0]

            if len( rt.products ) == 1:
                ### why this check, and why *2 here?
                if species1.radius * 2 < rt.products[0].radius:
                    raise RuntimeError,\
                        'radius of product must be smaller ' \
                        + 'than radius of reactant.'
            elif len( rt.products ) == 2:
                if species1.radius < rt.products[0].radius or\
                        species1.radius < rt.products[1].radius:
                    raise RuntimeError,\
                        'radii of both products must be smaller than ' \
                        + 'reactant.radius.'

            if self.reactionTypeMap1.has_key( species1 ):
                self.reactionTypeMap1[species1].append( rt )
            else:
                self.reactionTypeMap1[species1] = [ rt, ]

        elif numReactants == 2:
            species1 = rt.reactants[0]
            species2 = rt.reactants[1]
            self.reactionTypeMap2[ (species1,species2) ] = rt
            if species1 != species2:
                self.reactionTypeMap2[ (species2,species1) ] = rt

        else:
            raise RuntimeError, 'Invalid ReactionType.'


    def addInteractionType( self, it ):
        species = it.reactants[0]
        surface = it.reactants[1]
        self.interactionTypeMap[ (species, surface) ] = it


    def setAllRepulsive( self ):
        for species1 in self.speciesList.values():
            for species2 in self.speciesList.values():
                try:
                    _ = self.reactionTypeMap2[ ( species1, species2 ) ]
                except:
                    self.reactionTypeMap2[ ( species1, species2 ) ] =\
                                            RepulsionReactionType( species1,\
                                                                   species2 )
        for species in self.speciesList.values():
            for surface in self.surfaceList:
                try:
                    _ = self.interactionTypeMap[ ( species, surface ) ]
                except:
                    self.interactionTypeMap[ ( species, surface ) ] =\
                            SurfaceRepulsionInteractionType( species,surface )
        

    def throwInParticles( self, species, n, surface ):
        log.info( 'throwing in %s %s particles' % ( n, species.id ) )

        for i in range( int( n ) ):
            while 1:
                position = surface.randomPosition()
                if self.checkOverlap( position, species.radius ):
                    break
                else:
                    log.info( '%d-th particle rejected.' %i )
            
            self.createParticle( species, position, surface )


    def placeParticle( self, species, pos, surface=None ):
        log.info( 'place %s particle at %s' % ( species.id, pos ) )
        pos = numpy.array( pos )
        radius = species.radius

        if not self.checkOverlap( pos, radius ):
            raise NoSpace, 'overlap check failed'

        particle = self.createParticle( species, pos, surface )
        return particle


    def createParticle( self, species, pos, surface ):
        newserial = species.newParticle( pos )

        if isinstance( surface, CuboidalSurface ) or surface == None:
            # Particle in cytoplasm.
            # For backward compatibility.
            surface = self.defaultSurface

        newparticle = Particle( species, serial=newserial, surface=surface )
        self.addToParticleMatrix( newparticle, pos )
        return newparticle

    def removeParticle( self, particle ):
        particle.species.pool.removeBySerial( particle.serial )
        self.removeFromParticleMatrix( particle )

    def moveParticle( self, particle, newpos ):
        particle.pos = newpos
        self.updateOnParticleMatrix( particle, newpos )


    def addToParticleMatrix( self, particle, pos ):
        self.particleMatrix.add( particle,
                                 pos, particle.species.radius )

    def removeFromParticleMatrix( self, particle ):
        self.particleMatrix.remove( particle )

    def updateOnParticleMatrix( self, particle, pos ):
        self.particleMatrix.update( particle,
                                    pos, particle.species.radius )


    # checkOverlap return true if there are no particles within the
    # particles radius.
    # In subSpaceSimulator called from: fireSingleReaction. Asserts in:
    # propagateSingle, firePair, burstSingle, breakUpPair.
    def checkOverlap( self, pos, radius, ignore=[] ):
        
        particles, _ = \
            self.particleMatrix.getNeighborsWithinRadiusNoSort( pos, radius )
        if len( particles ) == 0:
            return True

        if [ p for p in particles if p not in ignore ]:
            #log.error( "Overlap with particle %s at position %s" % (p, pos) )
            return False
        else:
            return True


        
    def getParticlesWithinRadius( self, pos, radius, ignore=[] ):
        particles, _ =\
            self.particleMatrix.getNeighborsWithinRadius( pos, radius )

        return [ p for p in particles if p not in ignore ]


    def getParticlesWithinRadiusNoSort( self, pos, radius, ignore=[] ): 
        particles, _ =\
            self.particleMatrix.getNeighborsWithinRadiusNoSort( pos, radius )

        return [ p for p in particles if p not in ignore ]


    def clear( self ):

        self.dtMax = self.dtLimit
        self.dt = self.dtLimit


        
    '''
    Get closest n Particles.

    When the optional argument speciesList is given, only Particles of
    species in the list are considered.  When speciesList is not given
    or is None, all species in the simulator are considered.
    
    This method returns a tuple ( neighbors, distances ), where neighbors
    is a list of Particle objects.
    '''

    # There is a problem with these 4 methods: Particles can not be indexed.  
    # And they are not used anyway. Maybe before with gfrd?
    # In all these methods we try to get a reference to to the original 
    # particle, since we don't just want the values that the particleMatrix 
    # returns?

    """
    def getNeighborParticles( self, pos, n=None ):
        n, d = self.particleMatrix.getNeighbors( pos, n )
        # This doesn't seem right.
        neighbors = [ Particle( i[0], i[1] ) for i in n ]
        return neighbors, d

    def getNeighborParticlesNoSort( self, pos, n=None ):
        ### getNeighborsNoSort() does not have parameter n.
        n, d = self.particleMatrix.getNeighborsNoSort( pos, n )
        neighbors = [ Particle( i[0], i[1] ) for i in n ]
        return neighbors, d


    def getClosestParticle( self, pos, ignore=[] ):

        neighbors, distances =\
            self.getNeighborParticles( pos, len( ignore ) + 1 )

        for i in range( len( neighbors ) ): 
            if neighbors[i] not in ignore:
                closest, distance = neighbors[i], distances[i]

                #assert not closest in ignore
                return closest, distance

        # default case: none left.
        return None, INF
        #return DummyParticle(), INF


    def checkSurfaces( self, speciesIndex1, particleIndex ):

        speciesList = self.speciesList.values()

        species = speciesList[ speciesIndex1 ]
        pos = species.pool.positions[ particleIndex ].copy()

        dist = [ surface.distance( pos ) for surface in self.surfaceList ]

        if len( dist ) == 0:
            return -1, 0.0

        idx = numpy.argmin( dist )
        dist = dist[idx]

        dt = ( dist - species.radius ) ** 2 / \
                ( self.H * self.H * 6.0 * species.D ) 
        
        return dt, idx
        
    """

    def checkParticleMatrix( self ):

        if self.worldSize != self.particleMatrix.worldSize:
            raise RuntimeError,\
                'self.worldSize != self.particleMatrix.worldSize'


        ### why use numpy.array() here?
        total = numpy.array( [ species.pool.size for species 
                               in self.speciesList.values() ] ).sum()


        ### Check number of particles.
        if total != self.particleMatrix.size:
            raise RuntimeError,\
                'total number of particles %d != self.particleMatrix.size %d'\
                % ( total, self.particleMatrix.size )

        ### Check positions and radiuses of particles.
        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                particle = Particle( species, index=i )
                pos, radius = self.particleMatrix.get( particle )

                if ( particle.pos - pos ).sum() != 0:
                    raise RuntimeError,\
                        'particleMatrix positions consistency broken'
                if particle.species.radius != radius:
                    raise RuntimeError,\
                        'particleMatrix radii consistency broken'

    def check( self ):

        self.checkParticleMatrix()


    def dumpPopulation( self ):
        buf = ''
        for species in self.speciesList.values():
            buf += species.id + ':' + str( species.pool.size ) + '\t'

        return buf
