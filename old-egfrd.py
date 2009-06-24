def getNeighbors( self, pos, radius=INF, ignore=[] ):

    shells, dists = self.sphereMatrix.getNeighbors( pos )

    seen = dict.fromkeys( ignore )
    neighbors = []
    distances = []

    #print shells, dists
    for i, shell in enumerate( shells ):
        if not shell[0] in seen:
            seen[ shell[0] ] = None
            neighbors.append(shell[0])
            distances.append(dists[i])
            if dists[i] > radius:
                return neighbors, distances

    return neighbors + [DummySingle()], numpy.concatenate( [ distances,
                                                             [INF] ] )



    def getClosestShell( self, pos, ignore=[] ):
        neighbors, distances = self.getNeighborShells( pos )

        for i, neighbor in enumerate( neighbors ):
            if neighbor not in ignore:
                return neighbor, distances[i]

        return None, INF

    
    def getClosestNObjs( self, pos, n=1, ignore=[] ):

        neighbors, distances = self.getNeighborShells( pos, len( ignore ) + n )

        objs = []
        dists = []

        for i, neighbor in enumerate( neighbors ):
            if neighbor[0] not in ignore:
                objs += [neighbor[0]]
                dists += [distances[i]]
                if len( objs ) >= n:
                    return objs, dists

        return objs, dists


    def getNeighborShellsNoSort( self, pos, n=None ):
        return self.sphereMatrix.getNeighborsNoSort( pos, n )


    # Within radius, keys not extracted.
    def getClosestNObjs( self, pos, n=1, ignore=[] ):
    def getNeighborShellsWithinRadius( self, pos, radius ):
        return self.sphereMatrix.getNeighborsWithinRadius( pos, radius )


    def getNeighborShellsWithinRadiusNoSort( self, pos, radius ):
        return self.sphereMatrix.getNeighborsWithinRadiusNoSort( pos, radius )


    # With radius and ignore list.
    def getNeighborsWithinRadius( self, pos, radius, ignore=[] ):
        shells, distances =\
            self.sphereMatrix.getNeighborsWithinRadius( pos, radius )

        neighbors = [ s[0] for s in shells if s[0] not in ignore ]
        neighbors = uniq( neighbors )
        return neighbors
