import numpy

class Shape( object ):
    def distance( self, pos ):
        return abs( self.signedDistance( pos ) )


class Sphere( Shape ):
    def __init__( self, pos, size ):
        self.setParams( pos, size )

    def signedDistance( self, pos ):
        return math.sqrt( ( ( pos - self.pos ) ** 2 ).sum() ) - self.size

    def setParams( self, pos, size ):
        self.pos = numpy.array( pos )
        self.size = size

    def getParams( self ):
        return self.params

class Cylinder( Shape ):
    def __init__( self, pos, radius, orientation, size ):
        self.setParams( pos, radius, orientation, size )

    # Todo: this still assumes cylinder is infinitely long. Fix.
    def signedDistance( self, pos ):
        diff = numpy.linalg.norm( pos - self.projection(pos) )
        assert diff > self.radius
        distance = diff - self.radius
        print 'distance = ', distance
        return distance

    def setParams( self, pos, radius, orientation, size ):
        self.pos = numpy.array( pos )
        self.radius = radius
        # Todo: make sure orientation is normalized.
        self.orientation = numpy.array( orientation )
        self.size = numpy.array( size )

    def getParams( self ):
        return self.params

    # Return projection of 'pos' on the main axis of the cylinder.
    def projection( self, pos ):
        return self.pos + numpy.dot( (pos - self.pos), self.orientation ) * self.orientation
