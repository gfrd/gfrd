from singles import Single

# Todo.
class InteractionSingle( Single ):
    pass

# Todo.
class CylindricalSingle( InteractionSingle ):
    def __init__( self, pos, orientation, radius, halfLength ):
        self.pos = pos
        self.orientation = orientation
        self.radius = radius
        self.halfLength = halfLength

# Origin of box is where particles starts. Setting it at the corners is 1. not 
# easier, 2. mobilityRadius makes domain extent from (0+mobilityRadius to 
# L-mobilityRadius)
class CuboidalSingle( InteractionSingle ):
    def __init__( self, particle, reactionTypes, origin, baseVectors, widths, interactionRates ):
        self.origin = origin
        self.baseVectors = baseVectors
        self.widths = widths
        pos = toInternal( particle.pos )
        gf = OneDimGreensFunction( particle.species.D ) 
        Single.__init__( particle, reactionTypes, [ SimpleDomain(x, w, ks, gf) for (x, w, ks) in zip(pos, widths, interactionRates)] )

    def toInternal( self, pos ):
        # Todo.
        pass

    def toExternal( self, pos ):
        # Todo.
        pass

