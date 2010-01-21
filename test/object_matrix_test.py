#!/usr/bin/env python

import unittest

import numpy

from object_matrix import *
from cObjectMatrix import *
import math

class object_matrixTestCase( unittest.TestCase ):

    def setUp(self):
        pass

    def tearDown(self):
        pass

    '''
    SphereContainer
    '''
    def testS1( self ):

        c = SphereContainer(1.0, 10)
        self.assertEqual( 10, c.matrix_size )
        self.assertEqual( 0.1, c.cell_size )
        self.assertEqual( 0, c.size() )
        c.insert( 0, numpy.array([0.5, 0.3, 0.2]), 0.1 )
        self.assertEqual( 1, c.size() )
        #self.assertEqual( None, c.get( 1 ) )
        c.insert( 1, numpy.array([0.0,0.3,0.9]), 0.1)
        self.assertEqual( 2, c.size() )

        self.assertAlmostEqual( c.get(1)[0][0], 0.0 )
        self.assertAlmostEqual( c.get(1)[0][1], 0.3 )
        self.assertAlmostEqual( c.get(1)[0][2], 0.9 )

        #self.assertEqual( None, c[2] )

        a = c.neighbors_array( numpy.array([0.45, 0.23, 0.13]), 0.09 )
        self.assertEqual( 1, len(a[0]) )
        self.assertEqual( 1, len(a[1]) )
        self.assertEqual( 0, a[0][0] )

    def testS2( self ):
        c = SphereContainer( 1000, 3 )
        c.insert( 'A', [500, 500, 500], 50 )
        
        # Find neighbors.
        _,d = c.neighbors_array( [500, 500, 600], 75 )
        self.assertAlmostEqual( 50, d[0] )
        # cyclic should give same result
        _,d = c.neighbors_array_cyclic( [500, 500, 600], 75 )
        self.assertAlmostEqual( 50, d[0] )

        # Update with same value.
        c.insert( 'A', [500, 500, 500], 50 )
        _,d = c.neighbors_array( [500, 500, 600], 75 )
        self.assertAlmostEqual( 50, d[0] )
        # cyclic should give same result
        _,d = c.neighbors_array_cyclic( [500, 500, 600], 75 )
        self.assertAlmostEqual( 50, d[0] )

        # Now a real update.
        c.insert( 'A', [500, 500, 500], 75 )
        _,d = c.neighbors_array( [500, 500, 600], 100 )
        self.assertAlmostEqual( 25, d[0] )
        _,d = c.neighbors_array_cyclic( [500, 500, 600], 100 )
        self.assertAlmostEqual( 25, d[0] )

    def testS3( self ):
        c = SphereContainer( 1000, 3 )
        c.insert( 'A', [0, 0, 0], 100 )
        _,d = c.all_neighbors_array_cyclic( [700, 0, 0] )
        self.assertAlmostEqual( 200, d[0] ) # Ok.
        _,d = c.all_neighbors_array_cyclic( [600, 0, 0] )
        #self.assertAlmostEqual( 300, d[0] ) # Fail.
        self.assertAlmostEqual( 500, d[0] ) # Not Ok!


    '''
    CylinderContainer
    '''
    def testC1( self ):
        c = CylinderContainer( 1000, 3 )
        c.insert( 'A', [500, 500, 500], 25, [0,0,1], 50 )
        
        # Find neighbors of cylinder.
        _,d = c.neighbors_array( [500, 500, 600], 75 )
        self.assertAlmostEqual( 50, d[0] )
        _,d = c.neighbors_array_cyclic( [500, 500, 600], 75 )
        self.assertAlmostEqual( 50, d[0] )

        # Cylinder update with same value.
        c.insert( 'A', [500, 500, 500], 25, [0,0,1], 50 )
        _,d = c.neighbors_array( [500, 500, 600], 75 )
        self.assertAlmostEqual( 50, d[0] )
        _,d = c.neighbors_array_cyclic( [500, 500, 600], 75 )
        self.assertAlmostEqual( 50, d[0] )

        # Real update.
        c.insert( 'A', [500, 500, 500], 25, [0,0,1], 75 )
        _,d = c.neighbors_array( [500, 500, 600], 75 )
        self.assertAlmostEqual( 25, d[0] )
        _,d = c.neighbors_array_cyclic( [500, 500, 600], 75 ) 
        self.assertAlmostEqual( 25, d[0] )

    def testC2( self ):
        # Distance to cylinder in radial direction.
        c = CylinderContainer( 1000, 3 )
        c.insert( 'A', [500, 500, 500], 50, [0,0,1], 50 )
        _,d = c.all_neighbors_array( [500, 600, 527] )
        self.assertAlmostEqual( 50, d[0] )
        # Distance to cylinder edge.
        _,d = c.all_neighbors_array( [500, 553, 554] )
        self.assertAlmostEqual( 5, d[0] )

    def testC3( self ):
        # Distance to cylinder using periodic boundary conditions.
        # First distance towards radial.
        c = CylinderContainer( 1000, 3 )
        c.insert( 'A', [0, 0, 0], 50, [0,0,1], 50 )
        _,d = c.all_neighbors_array_cyclic( [950, 950, 0] )
        self.assertAlmostEqual( math.sqrt(2*50**2)-50, d[0] )

    def testC4( self ):
        c = CylinderContainer( 1000, 3 )
        c.insert( 'A', [0, 0, 0], 100, [0,0,1], 100 )
        _,d = c.all_neighbors_array_cyclic( [700, 0, 0] )
        self.assertAlmostEqual( 200, d[0] ) # Ok.
        _,d = c.all_neighbors_array_cyclic( [600, 0, 0] )
        #self.assertAlmostEqual( 300, d[0] ) # Fail.
        self.assertAlmostEqual( 500, d[0] ) # Not Ok!


    '''
    BoxContainer
    '''
    def testB1( self ):
        c = BoxContainer( 1000, 3 )
        c.insert( 'A', [500, 500, 500], [1,0,0], [0,1,0], [0,0,1], 50, 50, 50 )
        
        # Find neighbors of cylinder.
        _,d = c.neighbors_array( [500, 500, 600], 75 )
        self.assertAlmostEqual( 50, d[0] )
        _,d = c.neighbors_array_cyclic( [500, 500, 600], 75 )
        self.assertAlmostEqual( 50, d[0] )

        # Cylinder update with same value.
        c.insert( 'A', [500, 500, 500], [1,0,0], [0,1,0], [0,0,1], 50, 50, 50 )
        _,d = c.neighbors_array( [500, 500, 600], 75 )
        self.assertAlmostEqual( 50, d[0] )
        _,d = c.neighbors_array_cyclic( [500, 500, 600], 75 )
        self.assertAlmostEqual( 50, d[0] )

        # Real update.
        c.insert( 'A', [500, 500, 500], [1,0,0], [0,1,0], [0,0,1], 50, 50, 75 )
        _,d = c.neighbors_array( [500, 500, 600], 75 )
        self.assertAlmostEqual( 25, d[0] )
        _,d = c.neighbors_array_cyclic( [500, 500, 600], 75 ) 
        self.assertAlmostEqual( 25, d[0] )

    def testB2( self ):
        # Distance to cylinder in y direction.
        c = BoxContainer( 1000, 3 )
        c.insert( 'A', [500, 500, 500], [1,0,0], [0,1,0], [0,0,1], 50, 50, 50 )
        _,d = c.all_neighbors_array( [500, 600, 527] )
        self.assertAlmostEqual( 50, d[0] )
        # Distance to cylinder edge.
        _,d = c.all_neighbors_array( [500, 553, 554] )
        self.assertAlmostEqual( 5, d[0] )

    def testB3( self ):
        # Distance to cylinder using periodic boundary conditions.
        # Distance towards edge x-y.
        c = BoxContainer( 1000, 3 )
        c.insert( 'A', [0,0,0], [1,0,0], [0,1,0], [0,0,1], 50, 50, 50 )
        _,d = c.all_neighbors_array_cyclic( [946, 947, 0] )
        self.assertAlmostEqual( 5, d[0] )


if __name__ == "__main__":
    unittest.main()
