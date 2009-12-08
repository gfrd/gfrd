from surface import *

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


