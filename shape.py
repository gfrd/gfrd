import numpy

debug = 1

class Shape( object ):
    def __init__( self, distFunc ):
        # Note: we can use the distance function to compute the distance 
        # between 2 points, taken into account the periodic boundary 
        # conditions. Very cool. 
        self.distance = distFunc


    # DEBUG:
    def setPos ( self, pos ):
        raise RuntimeError, 'Shell doesnt have a pos anymore'
    def getPos ( self ):
        raise RuntimeError, 'Shell doesnt have a pos anymore'
    pos = property( getPos, setPos )


    def distanceTo( self, pos ):
        return abs( self.signedDistanceTo( pos ) )



class Sphere( Shape ):
    def __init__( self, origin, size, distFunc=None ):
        Shape.__init__( self, distFunc ) 
        # Todo: before we had a pos.copy() here. Think about it.
        self.origin = numpy.array( origin )
        self.size = size


    def signedDistanceTo( self, pos ):
        return self.distance( pos, self.origin ) - self.size



class Cylinder( Shape ):
    def __init__( self, origin, radius, orientation, size, distFunc=None ):
        Shape.__init__( self, distFunc )
        # Pos is the centre of the cylinder! This has the slight advantage 
        # over other options that we can make use of symmetry sometimes.
        self.origin = numpy.array( origin )
        self.radius = radius
        self.orientation = numpy.array( orientation ) / numpy.linalg.norm( orientation )
        # Size is the half length of the cylinder!
        self.size = numpy.array( size ) 


    # Look mom: no caps!
    def signedDistanceTo( self, pos ):
        posVector, zUnitVector, z, rUnitVector, r = self.toInternal( pos )
        zAbs = abs(z)

        if zAbs > self.size:
            # pos is (either) to the right or to the left of the cylinder.
            if r < self.radius:
                distance = zAbs
            else:
                # Difficult case. Compute distance to edge.
                edgVector = self.size * zUnitVector + self.radius * rUnitVector
                distance = numpy.norm( posVector - edgeVector )
        else:
            # pos is somewhere 'parellel' to the cylinder.
            distance = r - self.radius
        assert distance > 0
        return distance


    # Return the (z,r) components of pos in a coordinate system defined by the 
    # vectors zVector and rVector, where zVector is the orientation of the 
    # cylinder and rVector is choosen such that zVector and rVector define a 
    # plane in which pos lies.
    def toInternal( self, pos ):
        posVector = pos - self.origin

        z = numpy.dot( posVector, self.orientation )
        zUnitVector = self.orientation # By definition.
        zVector = z*zUnitVector

        rVector = posVector - zVector
        r = numpy.norm( rVector )
        rUnitVector = rVector / r

        assert pos == self.origin + zVector + rVector
        return posVector, zUnitVector, z, rUnitVector, r


    def __str__( self ):
        return "Cylinder: " + str( self.origin ) + " " + str( self.radius ) + " " + \
                              str( self.orientation ) + " " + str( self.size )
