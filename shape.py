import numpy
from utils import length, normalize

debug = 1

class Shape( object ):
    def __init__( self, distFunc ):
        # Note: we can use the distance function to compute the distance 
        # between 2 points, taken into account the periodic boundary 
        # conditions.
        self.distance = distFunc


    # DEBUG:
    def setPos ( self, pos ):
        raise RuntimeError, 'Shape doesnt have a pos anymore'
    def getPos ( self ):
        raise RuntimeError, 'Shape doesnt have a pos anymore'
    pos = property( getPos, setPos )


    def distanceTo( self, pos ):
        return abs( self.signedDistanceTo( pos ) )



class Sphere( Shape ):
    def __init__( self, origin, radius, distFunc=None ):
        Shape.__init__( self, distFunc ) 
        self.origin = numpy.array( origin )
        self.radius = radius


    def signedDistanceTo( self, pos ):
        return self.distance( pos, self.origin ) - self.radius



class Cylinder( Shape ):
    def __init__( self, origin, radius, orientationZ, size, distFunc=None ):
        Shape.__init__( self, distFunc )
        # Origin is the centre of the cylinder! This has the slight advantage 
        # over other options that we can make use of symmetry sometimes.
        self.origin = numpy.array( origin )
        self.radius = radius
        self.unitZ = normalize( numpy.array( orientationZ ) ) # Formerly known as orientationZ.
        # Size is the half length of the cylinder!
        self.size = size                                             # Formerly known as Lz.
        self.vectorZ = self.unitZ * size # Extra.

    def getSize( self ):
        return self._size
    def setSize( self, size ):
        self._size = size
    size = property( getSize, setSize )

    def signedDistanceTo( self, pos ):
        r, z, _, _ = self.toInternal( pos )
        dz = abs(z) - self.size
        dr = r - self.radius
        if dz > 0:
            # pos is (either) to the right or to the left of the cylinder.
            if dr > 0:
                # Compute distance to edge.
                distance = sqrt( dz*dz + dr*dr )
            else:
                distance = dz
        else:
            if dr > 0:
                # pos is somewhere 'parallel' to the cylinder.
                distance = dr
            else:
                # Inside cylinder, dz and dr are negative.
                distance = max( dz, dr )
        return distance


    # Return the (z,r) components of pos in a coordinate system defined by the 
    # vectors unitR and unitZ, where unitR is choosen such that unitR and 
    # unitZ define a plane in which pos lies.
    def toInternal( self, pos ):
        posVector = pos - self.origin

        z = numpy.dot( posVector, self.unitZ ) # Can be <0.
        posVectorZ = z*self.unitZ

        posVectorR = posVector - posVectorZ
        r = length( posVectorR )       # Always >= 0.
          
        return r, z, posVectorR, posVectorZ


    '''
    Returns:
    1. a vector pointing from the global origin to the projected point of 
       'pos' onto the main axis of the cylinder.
    2. a vector pointing from that point to 'pos'.
    3. the length of that 2nd vector.
    '''
    def calculateProjectionVectors( self, pos ):
        r, z, posVectorR, posVectorZ = self.toInternal( pos )
        return self.origin + posVectorZ, posVectorR, r


    def __str__( self ):
        return "Cylinder: " + str( self.origin ) + " " + str( self.radius ) + " " + \
                              str( self.unitZ ) + " " + str( self.size )


class DummyCylinder( Cylinder ):
    def __init__( self ):
        Cylinder.__init__( self, [0,0,0], 0, [1,1,1], 0 ) 


class Box( Shape ):
    def __init__( self, origin, vectorX, vectorY, vectorZ, Lx, Ly, Lz, distFunc=None ):
        Shape.__init__( self, distFunc )
        self.origin = numpy.array(origin)

        self.unitX = normalize(numpy.array(vectorX))
        self.unitY = normalize(numpy.array(vectorY))
        self.unitZ = normalize(numpy.array(vectorZ))

        assert Lx > 0 and Ly > 0 and Lz > 0
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz

        # vectorX, vectorY and vectorZ didn't have to be normalized before.
        self.vectorX = self.unitX * Lx
        self.vectorY = self.unitY * Ly
        self.vectorZ = self.unitZ * Lz


    def signedDistanceTo( self, pos ):
        x, y, z = self.toInternal( pos )
        dx = abs(x) - self.Lx
        dy = abs(y) - self.Ly
        dz = abs(z) - self.Lz

        if dx > 0:
            if dy > 0:
		if dz > 0:
		    distance = sqrt( dx*dx + dy*dy + dz*dz )
		else:
		    distance = sqrt( dx*dx + dy*dy )
	    else:
		if dz > 0:
		    distance = sqrt( dx*dx + dz*dz )
		else:
		    distance = dx
	else:
            if dy > 0:
		if dz > 0:
		    distance = sqrt( dy*dy + dz*dz )
		else:
		    distance = dy
	    else:
		if dz > 0:
		    distance = dz
		else:
                    # Inside box. Pick negative distance closest to 0.
		    distance = max( max(dx,dy), dz )

        return distance


    def toInternal( self, pos ):
        # First compute the (x,y,z) components of pos in a coordinate system 
        # defined by the vectors unitX, unitY, unitZ of the box.
        posVector = pos - self.origin
        x = numpy.dot( posVector, self.unitX )
        y = numpy.dot( posVector, self.unitY )
        z = numpy.dot( posVector, self.unitZ )

        return x, y, z


    '''
    Returns:
    1. a vector pointing from the global origin to the projected point of 
       'pos' onto the xy-plane of the box (or surface).
    2. a vector pointing from that point to 'pos'.
    3. the length of that 2nd vector.
    '''
    def calculateProjectionVectors( self, pos ):
        x, y, z = self.toInternal( pos )
        return self.origin + x * self.unitX + y * self.unitY, z * self.unitZ, z


    def __str__( self ):
        return "Box: " + str( self.origin ) + " " + str( self.vectorX ) + " " + \
                              str( self.vectorY ) + " " + str( self.vectorZ )

class DummyBox( Box ):
    def __init__( self ):
        Box.__init__( self, [0,0,0], [1,0,0], [0,1,0], [0,0,1], 1, 1, 1 ) 

