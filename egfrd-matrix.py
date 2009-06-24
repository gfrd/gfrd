


def setMatrixSize( self, size ):
    ParticleSimulatorBase.setMatrixSize( self, size )
    self.sphereMatrix.setMatrixSize( size )
    self.cylinderMatrix.setMatrixSize( size )


def getMatrixCellSize( self ):
    assert self.sphereMatrix.cellSize == self.cylinderMatrix.cellSize
    return self.sphereMatrix.cellSize


def addToShellMatrix( self, obj ):
    for i, shell in enumerate( obj.shellList ):
        key = (obj, i)
        if isinstance( shell, Sphere ):
            self.sphereMatrix.add( key, shell.origin, shell.size )
        elif isinstance( shell, Cylinder ):
            self.cylinderMatrix.add( key, shell )
        else: raise KeyError, 'Objecttype does not exit'


def removeFromShellMatrix( self, obj ):
    for i, shell in enumerate( obj.shellList ):
        key = (obj, i)
        # Check for type of shell -> remove key (is that the way to do 
        # it?)
        if isinstance( shell, Sphere ):
            self.sphereMatrix.remove( key )
        elif isinstance( shell, Cylinder ):
            self.cylinderMatrix.remove( key )
        else: raise KeyError, 'Objecttype does not exit'


def updateShellMatrix( self, obj ):
    #print 'updateMatrix'
    for i, shell in enumerate( obj.shellList ):
        key = (obj, i)
        #print 'update shell ', shell
        if isinstance( shell, Sphere ):
            self.sphereMatrix.update( key, shell.origin, shell.size )
        elif isinstance( shell, Cylinder ):
            self.cylinderMatrix.update( key, shell )
        else: raise KeyError, 'Objecttype does not exit'


'''
Find closest n shells.

This method returns a tuple ( neighbors, distances ).
'''

# Sort yes/no.
# Within radius yes/no.
# With ignore list yes/no.

# No radius (i.e. all), no ignore list, keys not extracted.


def getClosestObj( self, pos, ignore=[] ):
    """
    Don't be smart and do:
        keys, distances, _, _ = self.getNeighbors( pos, ignore )
        return keys[0], distance[0]
    since that would unneccesarily iterate one more time over all 
    objects in self.getNeighbors. 
    """
    closestSingle = DummySingle()
    closestDistance = INF

    for objectMatrix in [ self.sphereMatrix, self.cylinderMatrix ]:
        keys, distances = objectMatrix.getNeighbors( pos )
        """
        Don't be clever and think you can do:
            keys, distances = objectMatrix.getNeighbors( pos, 1 )
        Because then you are ignoring the ignore list.
        """
        for i, key in enumerate( keys ):
            single = key[0]
            if single not in ignore and distances[i] < closestDistance:
                closestSingle, closestDistance = single, distances[i]
                # Found yet a closer single. Break out of inner for loop 
                # and check other objectMatrices.
                break   

    return closestSingle, closestDistance


def getNeighborsWithinRadiusNoSort( self, pos, radius, ignore=[] ):
    keys, _ =\
        self.sphereMatrix.getNeighborsWithinRadiusNoSort( pos, radius )
    keys2, _ =\
        self.cylinderMatrix.getNeighborsWithinRadiusNoSort( pos, radius )

    neighbors = uniq( [ key[0] for key in keys if key[0] not in ignore ] )
    neighbors2 = uniq( [ key[0] for key in keys2 if key[0] not in ignore ] )

    neighbors.extend( neighbors2 )
    return neighbors


# Returns all singles and the distances towards their *shell* if that shell is 
# within a distance radius of pos, plus the closest single and the distance 
# towards it shell that is just outside of radius.
def getNeighbors( self, pos, radius=INF, ignore=[] ):
    # Diccionaries are more efficient for lookups.
    seen = dict.fromkeys( ignore )
    neighbors = []
    distances = []

    closestSingle = DummySingle()
    closestDistance = INF

    for objectMatrix in [ self.sphereMatrix, self.cylinderMatrix ]:
        keys, dists = objectMatrix.getNeighbors( pos )

        for i, key in enumerate( keys ):
            single = key[0]
            if not single in seen:
                # Since a single can in theory have more than 1 shell, and 
                # for each shell there is an entry in the objectMatrix, we 
                # signal here that we have already seen this single.
                # This is a bug: seen[ key ] = None
                #log.debug( 'dists[%d] = %g' % (i, dists[i]) )
                #log.debug( 'radius = %g', radius )
                seen[ single ] = None
                if dists[i] > radius:
                    # This is a single that has a shell that is more than 
                    # radius away from pos. If it is closer than the 
                    # closest such one we found so far: store it.
                    # Always break out of the inner for loop now and check 
                    # the other objectMatrices.
                    if  dists[i] < closestDistance:
                        closestSingle = single
                        closestDistance = dists[i]
                        break
                    else:
                        break # Just to be sure you get the point.
                else:
                    neighbors.append( single )
                    distances.append( dists[i] )

    #print 'binnen: ', neighbors, distances
    #print 'buiten: ', closestSingle, closestDistance
    return neighbors, distances, closestSingle, closestDistance
