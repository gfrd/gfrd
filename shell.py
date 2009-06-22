
# Still used for Pair.
class Shell( object ):
    def __init__( self, pos, size ):
        self.pos = pos.copy()
        self.size = size

    # DEBUG:
    def setRadius( self, radius ):
        raise RuntimeError, 'Shell doesnt have a radius anymore'
    def getRadius( self ):
        raise RuntimeError, 'Shell doesnt have a radius anymore'
    radius = property( getRadius, setRadius )
