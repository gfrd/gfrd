#!/usr/bin/env python


import numpy
# Todo: I want to do something like this:
#from gfrdbase import log

### object_matrix is a module (library) that is written in C++ and
### pythonified using Boost library.
### Source: object_matrix/src/object_container.hpp
### Pythonify: object_matrix/boost.python/peer/ObjectContainer.hpp
###   and via object_matrix/boost.python/object_matrix_module.hpp
import object_matrix

from utils import *

debug=1

class ObjectMatrix( object ):

    def __init__( self ):

        self.worldSize = 1.0
        self.setMatrixSize( 3 )
        self.initialize()


    def setWorldSize( self, size ):
        self.worldSize = size
        self.initialize()


    def setMatrixSize( self, size ):

        if size < 3:
            raise RuntimeError,\
                'Size of distance cell matrix must be at least 3'
        self.matrixSize = size
        self.initialize()


    def getSize( self ):

        return self.impl.size()

    size = property( getSize )


    def getCellSize( self ):
        return self.impl.cell_size

    cellSize = property( getCellSize )


    def initialize( self ):

        self.impl = object_matrix.SphereContainer( self.worldSize, 
                                                   self.matrixSize )
        self.cylinderContainer = object_matrix.CylinderContainer( self.worldSize,
                                                   self.matrixSize )



    def clear( self ):
        self.initialize()


    ### Adds a single/pair to matrix.
    def add( self, key, pos, radius ):

        if( debug ):
            try:
                # For pretty printing. Singles/Pairs/Cylinders:
                object = key[0]
            except TypeError:
                # Particles:
                object = key

        assert radius < self.cellSize * .5
        assert not self.impl.contains( key )

        ### Calls insert() in ObjectContainer.hpp. 
        self.impl.insert( key, pos, radius )

    def addCylinder( self, key, cylinder ):

        if( debug ):
            try:
                object = key[0]
            except TypeError:
                object = key

        # Todo. Do something like this:
        # assert radius < self.cellSize * .5
        assert not self.cylinderContainer.contains( key )

        self.cylinderContainer.insert( key, cylinder.pos, cylinder.orientation, cylinder.radius, cylinder.halfLength )

    def remove( self, key ):

        assert self.impl.contains( key )

        self.impl.erase( key )


    def removeCylinder( self, cylinder ):
        key = (cylinder, 0)
        assert self.cylinderContainer.contains( key )

        self.cylinderContainer.erase( key )


    def update( self, key, pos, radius ):

        if( debug ):
            try:
                object = key[0]
            except TypeError:
                object = key

        assert radius < self.cellSize * .5
        assert self.impl.contains( key )

        self.impl.update( key, pos, radius )


    def get( self, key ):

        return self.impl.get( key )


    def getNeighborsCyclicNoSort( self, pos ):

        return self.impl.all_neighbors_array_cyclic( pos )


    ### Called via some steps from subSpaceSimulator.getNeighborsWithinRadiusNoSort().
    ### Returns neighbors and distances.
    def getNeighborsWithinRadiusNoSort( self, pos, radius ):

        assert radius < self.cellSize * .5
        return self.impl.neighbors_array_cyclic( pos, radius )


    ### Called via some steps from subSpaceSimulator.getNeighbors().
    def getNeighborsCyclic( self, pos, n=None ):

        neighbors, distances = self.impl.all_neighbors_array_cyclic( pos )
        topargs = distances.argsort()[:n]

        return neighbors.take( topargs ), distances.take( topargs )


    def getNeighborsWithinRadius( self, pos, radius ):

        assert radius < self.cellSize * .5

        neighbors, distances = \
            self.impl.neighbors_array_cyclic( pos, radius )

        topargs = distances.argsort()

        return neighbors.take( topargs ), distances.take( topargs )


    ### Called from subSpaceSimulator.getNeighbors().
    def getNeighbors( self, pos, n=None ):

        return self.getNeighborsCyclic( pos, n )


    def getNeighborsNoSort( self, pos ):

        return self.getNeighborsCyclicNoSort( pos )


    def check( self ):

        pass

