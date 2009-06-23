from shape import Sphere

# Still used for Pair.
class Shell( Sphere ):
    def __init__( self, pos, size ):
        Sphere.__init__( self, pos, size )


    # DEBUG:
    def setPos ( self, pos ):
        raise RuntimeError, 'Shell doesnt have a pos anymore'
    def getPos ( self ):
        raise RuntimeError, 'Shell doesnt have a pos anymore'
    pos = property( getPos, setPos )


    # DEBUG:
    def setRadius( self, radius ):
        raise RuntimeError, 'Shell doesnt have a radius anymore'
    def getRadius( self ):
        raise RuntimeError, 'Shell doesnt have a radius anymore'
    radius = property( getRadius, setRadius )
