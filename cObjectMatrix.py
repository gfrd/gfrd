#!/usr/bin/env python
import numpy
import object_matrix
from utils import *
from shape import Sphere, Cylinder
# Todo: I would like to do something like this.
#from gfrdbase import log


'''
object_matrix is a module (library) that is written in C++ and
pythonified using Boost library.

Source:
    object_matrix/src/object_container.hpp
Distance to shell functionality: 
    object_matrix/boost.python/peer/ObjectContainer.hpp
Pythonifier:
    object_matrix/boost.python/object_matrix_module.hpp
'''
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


    def getMatrixSize( self ):
        self.matrixSize


    def getSize( self ):
        return self.impl.size()
    size = property( getSize )


    def getCellSize( self ):
        return self.impl.cell_size
    cellSize = property( getCellSize )


    def clear( self ):
        self.initialize()


    def remove( self, key ):
        # Same for CylinderMatrix and SphereMatrix.
        assert self.impl.contains( key )
        self.impl.erase( key )


    '''
    Different getNeighbors methods. Return neighbors and distances.
    Options:
        * Sort yes/no.
        * Within radius?
            + yes: impl.neighbors_array_cyclic
            + no: impl.all_neighbors_array_cyclic
    These 4 methods return the distance from pos to the *shell*, not to the 
    origin.
    '''
    def getNeighbors( self, pos, n=None ):
        neighbors, distances = self.impl.all_neighbors_array_cyclic( pos )
        topargs = distances.argsort()[:n] # Actually works when n is None!
        return neighbors.take( topargs ), distances.take( topargs )


    def getNeighborsNoSort( self, pos ):
        return self.impl.all_neighbors_array_cyclic( pos )


    def getNeighborsWithinRadius( self, pos, radius ):
        assert radius < self.cellSize * .5

        neighbors, distances = \
            self.impl.neighbors_array_cyclic( pos, radius )

        topargs = distances.argsort()
        return neighbors.take( topargs ), distances.take( topargs )


    def getNeighborsWithinRadiusNoSort( self, pos, radius ):
        assert radius < self.cellSize * .5
        return self.impl.neighbors_array_cyclic( pos, radius )


    def check( self ):
        pass


    ''' 
    Get all objects from all cells in the matrix.
    Calling all_neighbors_array_cyclic once is not sufficient.
    Ugly! Hack! Todo! Fixme! Now!
    '''
    def getAll( self ):
        if not isinstance( self.worldSize, float ): 
            raise NotImplementedError
        cellSize = self.cellSize
        numSteps = self.matrixSize / 2
        stepSize = cellSize * 2 # Skip every other cell.
        offset = cellSize / 2   # Sample midpoint of cell.

        steps = [ offset + stepSize * i for i in range( numSteps ) ]
        points = [[x,y,z] for x in steps for y in steps for z in steps ] 

        seen = {} 
        for point in points:
            particles, _ = self.impl.all_neighbors_array_cyclic( point )
            for particle in particles:
                if particle not in seen:
                    seen[ particle ] = None

        return seen.keys()


# Used for particles as particleMatrix in ParticleSimulatorBase, and for 
# spherical shells in EGFRDSimulator.
class SphereMatrix( ObjectMatrix ):
    def __init__( self ):
        ObjectMatrix.__init__( self )


    def initialize( self ):
        self.impl = object_matrix.SphereContainer( self.worldSize, self.matrixSize )


    def add( self, key, pos, radius ):
        assert not self.impl.contains( key )
        self.impl.insert( key, pos, radius )


    def update( self, key, pos, radius ):
        assert self.impl.contains( key )
        self.impl.update( key, pos, radius )


    def get( self, key ):
        # There is not an equivalent for cylinders.
        return self.impl.get( key )


    def addOrUpdateValue( self, key, object ):
        assert radius < self.cellSize * .5
        self.impl.insert( key, object.pos, object.radius )


class CylinderMatrix( ObjectMatrix ):
    def __init__( self ):
        ObjectMatrix.__init__( self )


    def initialize( self ):
        self.impl = object_matrix.CylinderContainer( self.worldSize, self.matrixSize )


    def add( self, key, shell ):
        # Todo. Do something like this.
        #assert radius < self.cellSize * .5
        assert not self.impl.contains( key )
        self.impl.insert( key, shell.origin, shell.radius, shell.orientation, shell.size )


    def update( self, key, shell ):
        assert self.impl.contains( key )
        # object_matrix handles updates nicely.
        self.impl.add( key, shell.origin, shell.radius, shell.orientation, shell.size )


class BoxMatrix( ObjectMatrix ):
    def __init__( self ):
        ObjectMatrix.__init__( self )


    def initialize( self ):
        self.impl = object_matrix.BoxContainer( self.worldSize, self.matrixSize )


    def add( self, key, shell ):
        # Todo. Do something like this.
        #assert radius < self.cellSize * .5
        assert not self.impl.contains( key )
        self.impl.insert( key, shell.origin, shell.unitX, shell.unitY, shell.unitZ, shell.Lx, shell.Ly, shell.Lz )


    def update( self, key, shell ):
        assert self.impl.contains( key )
        # object_matrix handles updates nicely.
        self.impl.add( key, shell.origin, shell.unitX, shell.unitY, shell.unitZ, shell.Lx, shell.Ly, shell.Lz )



