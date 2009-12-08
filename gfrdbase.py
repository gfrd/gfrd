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

class NoSpace( Exception ):
    pass


class Species( object ):
    def __init__( self, id, D=None, radius=None, surface=None ):
        self.id = id
        self.D = D
        self.radius = radius
        self.surface = surface
        self.pool = ParticlePool() # Stores positions only.


    def newParticle( self, position ):
        serial = self.pool.newParticle( position )
        return serial


    def removeParticleByIndex( self, index ):
        self.pool.removeByIndex( index )


    def removeParticleBySerial( self, serial ):
        self.pool.removeBySerial( serial )


    def __str__( self ):
        return self.id


class DummySpecies( object ):
    """This is needed internally during initialization for the virtual product 
    of a decay or surface absorption reaction.

    """
    def __init__( self, surface ):
        self.radius = 0
        self.surface = surface


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
                raise ValueError( 'total product radii must be smaller than '
                                  'reactant radius * 2' )
        if self.k < 0:
            raise ValueError( 'k < 0' )


    def order( self ):
        return len( self.reactants )


    def __str__( self ):
        s = 'k=' + str( self.k ) + ': '

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


#Reaction Types.
class UnimolecularReactionType( ReactionType ):
    """A -> B. Deprecated.

    """
    def __init__( self, s1, p1, k ):
        ReactionType.__init__( self, [ s1, ], [ p1, ], k )


class DecayReactionType( ReactionType ):
    """A -> None. Deprecated.

    """
    def __init__( self, s1, k ):
        ReactionType.__init__( self, [ s1, ], [], k )


class BindingReactionType( ReactionType ):
    """A + B -> C. Deprecated.

    """
    def __init__( self, s1, s2, p1, k ):
        ReactionType.__init__( self, [ s1, s2 ], [ p1, ], k )
        # Todo. These were not used, were they?
        #D = s1.D + s2.D
        #sigma = s1.radius + s2.radius


class UnbindingReactionType( ReactionType ):
    """C -> A + B. Deprecated.

    """
    def __init__( self, s1, p1, p2, k ):
        ReactionType.__init__( self, [ s1, ], [ p1, p2 ], k )


class RepulsionReactionType( ReactionType ):
    """A + B is repulsive. Deprecated.

    All combinations of species for which no BindingReactionType is defined 
    are repulsive by default.

    """
    def __init__( self, s1, s2 ):
        ReactionType.__init__( self, [ s1, s2 ], [], 0.0 )
        # Todo. These were not used, were they?
        #D = s1.D + s2.D
        #sigma = s1.radius + s2.radius


class SurfaceBindingReactionType( ReactionType ):
    """A + Surface -> B_on_Surface. Deprecated.

    """
    def __init__( self, reactantSpecies, productSpecies,  k ):
        ReactionType.__init__( self, [ reactantSpecies ], 
                               [ productSpecies, ], k )


class SurfaceAbsorptionReactionType( ReactionType ):
    """A + Surface -> None. Deprecated.

    """
    pass


class SurfaceDirectBindingReactionType( ReactionType ):
    """A + B_on_Surface + Surface -> C_on_Surface. Deprecated.

    A + Surface should be repulsive.

    Todo. Not yet implemented.

    """
    def __init__( self, reactantSpecies1, reactantSpecies2, 
                  productSpecies,  k ):
        ReactionType.__init__( self, [ reactantSpecies1, reactantSpecies2 ],
                               [ productSpecies, ], k )


class SurfaceRepulsionReactionType( ReactionType ):
    """A + Surface is repulsive. Deprecated.

    When no SurfaceBindingReactionType is defined for a combination of 
    species and surface, they are repulsive by default.

    """
    def __init__( self, species, surface ):
        ReactionType.__init__( self, [ species, surface ], [], 0.0 )


class SurfaceUnbindingReactionType( ReactionType ):
    """A_on_Surface -> B. Deprecated.
    
    Surface unbinding. Poisson process.

    After unbinding from a surface the particle always ends up on the 
    defaultSurface (world) for now.

    """
    def __init__( self, reactantSpecies, productSpecies, k ):
        ReactionType.__init__( self, [ reactantSpecies ],
                               [ productSpecies, ], k )


class SurfaceDirectUnbindingReactionType( ReactionType ):
    """A_on_Surface -> B_on_Surface + C. Deprecated.

    After unbinding from a surface, particle2 always ends up on the 
    defaultSurface (world), and particle1 stays on the surface it was on.

    Todo. Not yet implemented.

    """
    def __init__( self, reactantSpecies, productSpecies1, productSpecies2, k ):
        ReactionType.__init__( self, [ reactantSpecies ], 
                               [ productSpecies1, productSpecies2 ], k )


class Reaction:
    def __init__( self, type, reactants, products ):
        self.type = type
        self.reactants = reactants
        self.products = products
    def __str__( self ):
        return ( 'Reaction( ' + str( self.type ) + ', ' + 
                 str( self.reactants ) + ', ' + str( self.products ) + ' )' )


class Particle( object ):
    def __init__( self, species, serial=None, index=None ):
        """Not only constructor, but can also be used to make a copy of a
        particle by index or serial.

        """
        self.species = species

        if not serial is None:
            self.serial = serial
        elif not index is None:
            self.serial = species.pool.getSerialByIndex( index )
        else:
            raise ValueError( 'give either serial or index.' )

        self.hash = hash( self.species ) ^ self.serial

        self.surface = species.surface


    def posString( self ):
        # Todo. Try catch block should not be needed.
        try: 
            pos = self.getPos()
            factor = 1
            return '(%.2g %.2g %.2g)' % \
                   ( pos[0] * factor, pos[1] * factor, pos[2] * factor ) 
        except KeyError:
            return 'None'


    def __str__( self ):
        return ( "( " + self.species.id  + ", " + str( self.serial ) +
                 ' ). pos = ' + self.posString() )


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
        return pool.positions[pool.indexMap[self.serial]]


    def setPos( self, newpos ):
        pool = self.species.pool
        pool.positions[pool.indexMap[self.serial]] = newpos
    pos = property( getPos, setPos )


    def getRadius( self ):
        return self.species.radius
    radius = property( getRadius )


    def getIndex( self ):
        return self.species.pool.indexMap[self.serial]


class DummyParticle( object ):
    def __init__( self ):
        self.species = None
        self.serial = -1



class ParticlePool( object ):
    """ParticlePool

    Stores positions only.

    """
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

        self.indexMap[newserial] = newindex
        self.serials[newindex] = newserial
        self.positions[newindex] = position

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

        serialOfRemovedItem = self.serials[index]
        serialOfMovedItem = self.serials[self.size]

        # swap items at index and the end, and discard the item
        # to be removed, which is now at the end, by shrinking
        # the arrays by one.
        self.serials[index] = self.serials[self.size]
        self.positions[index] = self.positions[self.size]
        self.__resizeArrays( self.size )

        # book keeping
        del self.indexMap[serialOfRemovedItem]
        self.indexMap[serialOfMovedItem] = index


    def getSerialByIndex( self, index ):
        return self.serials[index]


    def getIndex( self, serial ):
        return self.indexMap[serial]



class ParticleSimulatorBase( object ):
    def __init__( self, worldSize ):
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

        # Stores complete particles (unlike particlePool).
        self.particleMatrix = SphereMatrix()

        self.setWorldSize( worldSize )

        self.lastReaction = None

        # Particles of a Species whose surface is not specified will be add to 
        # the world.
        # The world has to be cubic, because the objectmatrix neeeds it to be.
        # This line here is also why worldSize has to be specified in the
        # constructor now, and should not be redefined.
        self.defaultSurface = \
            CuboidalRegion( [ 0, 0, 0 ], [ worldSize, worldSize, worldSize ], 
                             'world' )


    def initialize( self ):
        # Check if user specified a surface for each species that is a product 
        # of an interaction.
        '''Todo.
        for interaction in self.interactionTypeMap.itervalues():
            product = interaction.products[0]
            if product.surface == self.defaultSurface:
                    raise RuntimeError( 'Use addSpecies(' + product.id + 
                                        ', someSurface), because ' +
                                        product.id + ' seems to be ' +
                                        'the product of interaction: ' +
                                        str( interaction ) )
        '''


    def getClosestSurface( self, pos, ignore ):
        """Return sorted list of tuples with:
            - distance to surface
            - surface itself

        We can not use objectmatrix, it would miss a surface if the origin of the 
        surface would not be in the same or neighboring cells as pos.

        """
        surfaces = [ None ]
        distances = [ INF ]

        ignoreSurfaces = []
        for obj in ignore:
            if isinstance( obj.surface, Surface ):
                # For example ignore surface that particle is currently on.
                ignoreSurfaces.append( obj.surface )

        for surface in self.surfaceList:
            if surface not in ignoreSurfaces:
                posTransposed = cyclicTranspose( pos, surface.origin, 
                                                 self.worldSize )
                distanceToSurface = surface.signedDistanceTo( posTransposed )
                distances.append( distanceToSurface )
                surfaces.append( surface )

        return min( zip( distances, surfaces ) )


    def getClosestSurfaceWithinRadius( self, pos, radius, ignore ):
        """Return sorted list of tuples with:
            - distance to surface
            - surface itself

        """
        distanceToSurface, closestSurface = self.getClosestSurface( pos, 
                                                                    ignore ) 
        if distanceToSurface < radius:
            return distanceToSurface, closestSurface
        else:
            return INF, None


    def reconstructParticleMatrix( self ):
        self.particleMatrix.clear()
        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                particle = Particle( species, index=i )
                self.particleMatrix.add( particle, particle.pos, 
                                         species.radius )


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
            log.debug( '\t\tdistance_Simple.. is used' )
            self._distance = distance_Simple
            #self._distanceArray = distanceArray_Simple
            self._distanceSq = distanceSq_Simple
            self._distanceSqArray = distanceSqArray_Simple
        else:
            log.debug( '\t\tdistance_Cyclic.. is used' )
            self._distance = distance_Cyclic
            #self._distanceArray = distanceSqArray_Cyclic
            self._distanceSq = distanceSq_Cyclic
            self._distanceSqArray = distanceSqArray_Cyclic


    def getWorldSize( self ):
        return self.worldSize


    def setMatrixSize( self, size ):
        self.particleMatrix.setMatrixSize( size )


    def applyBoundary( self, pos ):
        """Use this method to account for periodic boundary conditions.
        
        """
        pos %= self.worldSize


    def getReactionType1( self, species ):
        return self.reactionTypeMap1.get( species, None )


    def getReactionType2( self, species1, species2 ):
        return self.reactionTypeMap2.get( ( species1, species2 ), None )


    def getInteractionType( self, species, surface ):
        return self.interactionTypeMap.get( ( species, surface ) )
    

    def getSpeciesByIndex( self, i ):
        return self.speciesList.values()[i]


    def distanceSq( self, position1, position2 ):
        return self._distanceSq( position1, position2, self.worldSize )


    def distance( self, position1, position2 ):
        return self._distance( position1, position2, self.worldSize )
        

    def distanceSqArray( self, position1, positions ):
        return self._distanceSqArray( position1, positions, self.worldSize )


    def distanceArray( self, position1, positions ):
        return numpy.sqrt( self.distanceSqArray( position1, positions ) )


    def addPlanarSurface( self, origin, vectorX, vectorY, Lx, Ly, Lz=0, 
                          name="PlanarSurface" ):
        """Add a planar surface.

        origin -- [ x0, y0, z0 ] is the *center* of the planar surface.
        vectorX -- [ x1, y1, z1 ] and
        vectorY -- [ x2, y2, z2 ] are 2 perpendicular vectors that don't have 
        to be normalized that span the plane. For example [1,0,0] and [0,1,0]
        for a plane at z=0.

        Lx -- lx and 
        Ly -- ly are the distances from the origin of the plane along vectorX 
            or vectorY *to an edge* of the plane. PlanarSurfaces are finite.
        Lz -- dz, the thickness of the planar surface, can be omitted for Lz=0.
        name -- a descriptive name, only for nicer output.

        """
        return self.addSurface( PlanarSurface( origin, vectorX, vectorY, Lx, 
                                               Ly, Lz, name ) )

    def addCylindricalSurface( self, origin, radius, orientation, size, 
                               name="CylindricalSurface" ):
        """Add a cylindrical surface.

        origin -- [ x0, y0, z0 ] is the *center* of the cylindrical surface.
        radius -- r is the radis of the cylinder.
        orientation -- [ x1, y1, z1 ] is a vector that doesn't have to
            normalized that defines the orienation of the cylinder. For 
            example [0,0,1] for a for a cylinder along the z-axis.
        size -- lz is the distances from the origin of the cylinder along 
            the oriention vector to the end of the cylinder. So effectively
            the *half-length*. CylindricalSurfaces are finite.
        name -- a descriptive name, only for nicer output.

        """
        return self.addSurface( CylindricalSurface( origin, radius, 
                                                    orientation, size, name ) )


    def addSurface( self, surface ):
        if ( not isinstance( surface, Surface ) or
             isinstance( surface, CuboidalRegion ) ):
            raise RuntimeError( str( surface ) + ' is not a surface.' )

        self.surfaceList.append( surface )
        return surface


    def addSpecies( self, species, surface=None, D=None, radius=None ):
        """Add a new species.

        A species is a type of particles. By default the species is added to 
        the 'world'. If a surface is specified, it is added to that surface. Per 
        surface a different diffusion constant D and radius can be specified. 
        By default the ones for the 'world' are used.  

        species -- a species created with Species().
        surface -- the surface this species can exist on.
        D       -- the diffusion constant of the particles.
        radius  -- the radii of the particles.

        """
        createNewSpecies = False

        if surface == None:
            surface = self.defaultSurface
        else:
            assert any( surface == s for s in self.surfaceList ), \
                   '%s not in surfaceList.' % ( surface )

            # Construct new id if species lives on a surface.
            id = species.id + str( surface )

            createNewSpecies = True

        if D == None:
            assert species.D != None, \
                   'Diffusion constant of species %s not specified: ' % species
            D = species.D
        else:
            createNewSpecies = True

        if radius == None:
            assert species.radius != None, \
                   'Radius of species not specified: ' % species
            radius = species.radius
        else:
            createNewSpecies = True

        if createNewSpecies:
            # Create a new species for internal use. Don't use the user 
            # defined species at all. The new id is a concatenation of the 
            # user defined id and the surface this species exists on.
            species = Species( id , D, radius, surface ) 
        else:
            # Don't create a new species, use the one the user had defined.  
            # Only set the surface to be the defaultSurface.
            species.surface = self.defaultSurface

        assert not self.speciesList.has_key( species.id ), \
               'Species with id = %s has already been added.' % ( species.id )

        self.speciesList[ species.id ] = species

        return species


    def isSurfaceBindingReaction( self, rt ):
        # A surface binding reaction that is a surface absorption reaction 
        # doesn't have a product species. This check always works, but don't 
        # call it before the interaction is added to the interactionTypeMap.
        return any( [ rt == it for it in self.interactionTypeMap.values() ] )


    def isSurfaceUnbindingReaction( self, rt ):
        currentSurface = rt.reactants[0].surface
        targetSurface = rt.products[0].surface

        return ( currentSurface != self.defaultSurface and
                 targetSurface == self.defaultSurface )


    def isDirectSurfaceBindingReaction( self, rt ):
        currentSurface = rt.reactants[0].surface
        targetSurface1 = rt.products[0].surface
        targetSurface2 = rt.products[1].surface

        return( currentSurface == self.defaultSurface and
                xor( targetSurface1 != self.defaultSurface,
                     targetSurface2 != self.defaultSurface ) )


    def isDirectSurfaceUnbindingReaction( self, rt ):
        currentSurface = rt.reactants[0].surface
        targetSurface1 = rt.products[0].surface
        targetSurface2 = rt.products[1].surface

        return( currentSurface != self.defaultSurface and
                xor( targetSurface1 == self.defaultSurface,
                     targetSurface2 == self.defaultSurface ) )


    def convertSpeciesSurfaceTupleToSpecies( self, tuple ):
        """Helper.

        """
        if isinstance( tuple, Species ):
            # So it wasn't a tuple.
            species = tuple
            surface = self.defaultSurface
            id = species.id
        else:
            # Unpack (species, surface)-tuple.
            species = tuple[0]
            surface = tuple[1]

            if species == 0:
                # This is the virtual product of a decay or surface absorption 
                # reaction.
                return DummySpecies( surface )

            # Note: see addSpecies for how id is constructed. 
            id = species.id + str( surface )

        try: 
            # Retrieve species from self.speciesList.
            return self.speciesList[ id ]
        except KeyError:
            raise RuntimeError( 'Species %s on surface %s '
                                'does not exist: ' % ( species, surface ) )


    def addReaction( self, reactants, products, k ): 
        reactants = map( self.convertSpeciesSurfaceTupleToSpecies, reactants )
        products  = map( self.convertSpeciesSurfaceTupleToSpecies, products )

        rt = ReactionType( reactants, products, k )

        if( len( rt.products ) == 1 ):
            currentSurface = rt.reactants[0].surface
            targetSurface = rt.products[0].surface

            if( currentSurface == self.defaultSurface and
                targetSurface != self.defaultSurface ):
                # Remove DummySpecies from products list. Were needed for 
                # surface absorption reaction.
                rt.products = [ product for product in rt.products if not
                                isinstance( product, DummySpecies ) ]
                # Surface binding is not a Poisson process, so this reaction 
                # should not be added to the reaction list but to the 
                # interaction list.
                self.interactionTypeMap[( rt.reactants[0], targetSurface )] = rt
                return

        # Remove DummySpecies from products list. Were needed for decay 
        # reaction.
        rt.products = [ product for product in rt.products if not
                        isinstance( product, DummySpecies ) ]
        self.addReactionType( rt )


    def addReactionType( self, rt ):
        numReactants = len( rt.reactants )

        if numReactants == 1:
            species1 = rt.reactants[0]

            if len( rt.products ) == 1:
                # Todo. Why this check, and why *2 here?
                if species1.radius * 2 < rt.products[0].radius:
                    raise RuntimeError( 'radius of product must be smaller '
                                        'than radius of reactant.' )
            elif len( rt.products ) == 2:
                if ( species1.radius < rt.products[0].radius or
                     species1.radius < rt.products[1].radius ):
                    raise RuntimeError( 'radii of both products must be '
                                        'smaller than reactant.radius.' )

            if self.reactionTypeMap1.has_key( species1 ):
                self.reactionTypeMap1[species1].append( rt )
            else:
                self.reactionTypeMap1[species1] = [ rt, ]

        elif numReactants == 2:
            species1 = rt.reactants[0]
            species2 = rt.reactants[1]
            self.reactionTypeMap2[( species1, species2 )] = rt
            if species1 != species2:
                self.reactionTypeMap2[( species2, species1 )] = rt

        else:
            raise RuntimeError( 'Invalid ReactionType.' )


    def setAllRepulsive( self ):
        for species1 in self.speciesList.values():
            for species2 in self.speciesList.values():
                try:
                    _ = self.reactionTypeMap2[( species1, species2) ]
                except:
                    self.reactionTypeMap2[( species1, species2 )] = \
                        ReactionType( [ species1, species2 ], [], 0.0 )

        for species in self.speciesList.values():
            for surface in self.surfaceList:
                try:
                    _ = self.interactionTypeMap[( species, surface )]
                except:
                    self.interactionTypeMap[( species, surface )] = \
                        ReactionType( [ species, surface ], [], 0.0 )


    def throwInParticles( self, species, n, surface=None ):
        if surface == None or isinstance( surface, CuboidalRegion ):
            # Species that is used internally is the same, but still some 
            # checks are done which are usefull.
            species = self.convertSpeciesSurfaceTupleToSpecies( species )  
        else:
            # Look up species that is used internally.
            species = self.convertSpeciesSurfaceTupleToSpecies( (species,
                                                                 surface ) )  

        surface = species.surface
        log.info( '\tthrowing in %s %s particles to %s' % ( n, species.id, 
                                                            surface ) )

        i = 0
        while i < int( n ):
            position = surface.randomPosition()

            # Check overlap.
            if self.checkOverlap( position, species.radius ):
                create = True
                # Check if not too close to a neighbouring surfaces for 
                # particles added to the world, or added to a self-defined 
                # box.
                if surface == self.defaultSurface or \
                   ( surface != self.defaultSurface and 
                     isinstance( surface, CuboidalRegion ) ):
                    distance, closestSurface = self.getClosestSurface( position,
                                                                       [] )
                    if ( closestSurface and
                         distance < closestSurface.minimalDistanceFromSurface( 
                                    species.radius ) ):
                        log.info( '\t%d-th particle rejected. To close to '
                                  'surface. I will keep trying.' % i )
                        create = False
                if create:
                    # All checks passed. Create particle.
                    self.createParticle( species, position )
                    i += 1
            else:
                log.info( '\t%d-th particle rejected. I will keep trying.' % i )


    def placeParticle( self, species, pos ):
        log.info( '\tplace %s particle at %s' % ( species.id, pos ) )
        pos = numpy.array( pos )
        radius = species.radius

        if not self.checkOverlap( pos, radius ):
            raise NoSpace( 'overlap check failed' )

        # Look up species that is used internally.
        species = self.convertSpeciesSurfaceTupleToSpecies( species )  

        particle = self.createParticle( species, pos )
        return particle


    def createParticle( self, species, pos ):
        newserial = species.newParticle( pos )

        newparticle = Particle( species, serial=newserial )
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


    def checkOverlap( self, pos, radius, ignore=[] ):
        """Return true if there are no particles within the particles radius.

        In subSpaceSimulator called from: fireSingleReaction. Asserts in:
        propagateSingle, firePair, burstSingle, breakUpPair.

        """
        particles, _ = \
            self.particleMatrix.getNeighborsWithinRadiusNoSort( pos, radius )
        if len( particles ) == 0:
            return True

        if [ p for p in particles if p not in ignore ]:
            #log.error( "Overlap with particle %s at position %s" % ( p, pos ) )
            return False
        else:
            return True

        
    def getParticlesWithinRadius( self, pos, radius, ignore=[] ):
        particles, _ = \
            self.particleMatrix.getNeighborsWithinRadius( pos, radius )

        return [ p for p in particles if p not in ignore ]


    def getParticlesWithinRadiusNoSort( self, pos, radius, ignore=[] ): 
        particles, _ = \
            self.particleMatrix.getNeighborsWithinRadiusNoSort( pos, radius )

        return [ p for p in particles if p not in ignore ]


    def clear( self ):
        self.dtMax = self.dtLimit
        self.dt = self.dtLimit


        

    # There is a problem with these 4 methods: Particles can not be indexed.  
    # And they are not used anyway. Maybe before with gfrd?
    # In all these methods we try to get a reference to to the original 
    # particle, since we don't just want the values that the particleMatrix 
    # returns?
    '''
    def getNeighborParticles( self, pos, n=None ):
        """
        Get closest n Particles.

        When the optional argument speciesList is given, only Particles of
        species in the list are considered.  When speciesList is not given
        or is None, all species in the simulator are considered.
        
        This method returns a tuple ( neighbors, distances ), where neighbors
        is a list of Particle objects.

        """
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
        neighbors, distances = self.getNeighborParticles( pos, 
                                                          len( ignore ) + 1 )

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

        species = speciesList[speciesIndex1]
        pos = species.pool.positions[particleIndex].copy()

        dist = [ surface.distance( pos ) for surface in self.surfaceList ]

        if len( dist ) == 0:
            return -1, 0.0

        idx = numpy.argmin( dist )
        dist = dist[idx]

        dt = ( dist - species.radius ) ** 2 / \
             ( self.H * self.H * 6.0 * species.D ) 
        
        return dt, idx
    '''


    def checkParticleMatrix( self ):
        if self.worldSize != self.particleMatrix.worldSize:
            raise RuntimeError( 'self.worldSize != '
                                'self.particleMatrix.worldSize' )


        ### why use numpy.array() here?
        total = numpy.array( [ species.pool.size for species 
                               in self.speciesList.values() ] ).sum()


        ### Check number of particles.
        if total != self.particleMatrix.size:
            raise RuntimeError( 'total number of particles %d != '
                                'self.particleMatrix.size %d' % 
                                ( total, self.particleMatrix.size ) )

        ### Check positions and radiuses of particles.
        for species in self.speciesList.values():
            for i in range( species.pool.size ):
                particle = Particle( species, index=i )
                pos, radius = self.particleMatrix.get( particle )

                if ( particle.pos - pos ).sum() != 0:
                    raise RuntimeError( 'particleMatrix positions consistency '
                                        'broken' )
                if particle.species.radius != radius:
                    raise RuntimeError( 'particleMatrix radii consistency '
                                        'broken' )


    def check( self ):
        self.checkParticleMatrix()


    def dumpPopulation( self ):
        buf = ''
        for species in self.speciesList.values():
            buf += species.id + ':' + str( species.pool.size ) + '\t'

        return buf


    def dumpReactions( self ):
        buf = 'Monomolecular reactions:\n'
        for reactions in self.reactionTypeMap1.itervalues():
            for reaction in reactions:
                buf += str( reaction ) + '\n'

        buf += '\nReactions of 2 particles:\n'
        # Select unique reactions.
        reactions = set( self.reactionTypeMap2.itervalues() )
        for reaction in reactions:
            # Ignore reflecting.
            if reaction.k > 0:
                buf += str( reaction ) + '\n'

        buf += '\nInteractions between a particle and a surface:\n'
        for (_, surface), interaction in self.interactionTypeMap.iteritems():
            # Ignore reflecting.
            if interaction.k > 0:
                buf += str( interaction ) + " (" + str( surface ) + ")" + '\n'
        return buf
