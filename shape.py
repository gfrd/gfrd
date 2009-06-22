import numpy

class Shape( object ):
    def distance( self, pos ):
        return abs( self.signedDistance( pos ) )


class Sphere( Shape ):
    def __init__( self, pos, size ):
        self.pos = numpy.array( pos )
        self.size = size

    def signedDistance( self, pos ):
        return math.sqrt( ( ( pos - self.pos ) ** 2 ).sum() ) - self.size


class Cylinder( Shape ):
    def __init__( self, pos, radius, orientation, size ):
        self.pos = numpy.array( pos )           # Middle!
        self.radius = radius
        # Todo: make sure orientation is normalized.
        self.orientation = numpy.array( orientation )
        self.size = numpy.array( size )         # half length

    # Todo: this still assumes cylinder is infinitely long. Fix.
    def signedDistance( self, pos ):
        diff = numpy.linalg.norm( pos - self.projection(pos) )
        distance = diff - self.radius
        assert distance > 0
        return distance

    # Return projection of 'pos' on the main axis of the cylinder.
    def projection( self, pos ):
        return self.pos + numpy.dot( (pos - self.pos), self.orientation ) * self.orientation

    def __str__( self ):
        return "Cylinder: " + str( self.pos ) + " " + str( self.radius ) + " " + \
                              str( self.orientation ) + " " + str( self.size )
