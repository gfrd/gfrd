#!/usr/bin/env python

__author__    = 'Laurens Bossen'
__license__   = 'Don\'t touch my shit'
__copyright__ = ''


import unittest

import _gfrd as mod

import numpy


class FirstPassageGreensFunction1DRadTestCase( unittest.TestCase ):

    def setUp( self ):
        pass

    def tearDown( self ):
        pass
    
    def test_Instantiation( self ):
        D = 1e-12
        kf = 1e-8
        a = 1e-7

        gf = mod.FirstPassageGreensFunction1DRad( D, kf )
        self.failIf( gf == None )
        gf.seta( a )


    def test_DrawTime( self ):
        D = 1e-12
        kf = 1e-8
        a = 1e-7
        r0 = 5e-8
        
        gf = mod.FirstPassageGreensFunction1DRad( D, kf  )
        gf.seta( a )

	
        t = gf.drawTime( 0.5, r0 )
        self.failIf( t <= 0.0 or t >= numpy.inf )

        t = gf.drawTime( 0.0, r0 )
        self.failIf( t < 0.0 or t >= numpy.inf )
	

        t = gf.drawTime( 1 - 1e-16, r0 )
        self.failIf( t <= 0.0 or t >= numpy.inf )


    def test_DrawTime_a_equal_sigma( self ):
        D = 1e-12
        kf = 1e-8
        a = 0
        r0 = a
        
        gf = mod.FirstPassageGreensFunction1DRad( D, kf  )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        self.assertEqual( 0.0, t )

    def test_DrawTime_a_near_sigma( self ):
        D = 1e-12
        kf = 1e-8
        a = 1e-14
        r0 = 0
        
        gf = mod.FirstPassageGreensFunction1DRad( D, kf  )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        self.failIf( t <= 0.0 or t >= numpy.inf )

    def test_DrawTime_r0_equal_a( self ):
        D = 1e-12
        kf = 1e-8
        a = 1e-7
        r0 = a
        
        gf = mod.FirstPassageGreensFunction1DRad( D, kf )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        self.assertEqual( 0.0, t )

    def test_DrawTime_r0_equal_sigma_kf_zero( self ):
        D = 1e-12
        kf = 0.0 # note this
        a = 1e-7
        r0 = -a
        
        gf = mod.FirstPassageGreensFunction1DRad( D, kf  )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        self.failIf( t < 0.0 or t >= numpy.inf )


    def no_test_DrawTime_r0_equal_sigma_kf_large( self ):
        D = 1e-12
        kf = 1e-8
        a = 10e-7
        r0 = -a + 1e-12
        
        gf = mod.FirstPassageGreensFunction1DRad( D, kf  )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        self.failIf( t < 0.0 or t >= numpy.inf )


    def test_DrawEventType( self ):
        D = 1e-12
        kf = 1e-8
        a = 1e-7
        r0 = 0

        gf = mod.FirstPassageGreensFunction1DRad( D, kf  )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        eventType = gf.drawEventType( 0.5, r0, t )
        self.failIf( eventType != 0 and eventType != 1 and eventType != 2 )

        eventType = gf.drawEventType( 0.0, r0, t )
        self.assertEqual( eventType, 0 )

        eventType = gf.drawEventType( 0.999999, r0, t )
        self.assertEqual( eventType, 1 )


    def no_test_DrawEventType_smallt( self ):
        D = 1e-12
        kf = 1e-8
        a = 1e-6 
        r0 = 0

        gf = mod.FirstPassageGreensFunction1DRad( D, kf  )
        gf.seta( a )

        t = gf.drawTime( 0.001, r0 )

        eventType = gf.drawEventType( 0.5, r0, t )
        self.failIf( eventType != 0 and eventType != 1 and eventType != 2 )

        eventType = gf.drawEventType( 0.0, r0, t )
        self.assertEqual( eventType, 0 )

        eventType = gf.drawEventType( 0.9999, r0, t )
        self.assertEqual( eventType, 1 )


    def test_DrawR( self ):
        D = 1e-12
        kf = 1e-8
        a = 1e-7
        r0 = 0
        
        gf = mod.FirstPassageGreensFunction1DRad( D, kf  )
        gf.seta( a )

        t = 1e-3

        r = gf.drawR( 0.5, r0, t )
        self.failIf( r < -a or r > a )

        r1 = gf.drawR( 0.0, r0, t )
        r2 = gf.drawR( 0.999999999999, r0, t )

        self.failIf( r1 < -a or r1 > a )
        self.failIf( r2 < -a or r2 > a )

        self.assertAlmostEqual( r1, -a )
        self.assertAlmostEqual( r2, a )


    def test_DrawR_zerot( self ):
        D = 1e-12
        kf = 1e-8
        a = 1e-7
        r0 = 0
        
        gf = mod.FirstPassageGreensFunction1DRad( D, kf  )
        gf.seta( a )

        t = 0.0

        r = gf.drawR( 0.5, r0, t )
        self.assertEqual( r0, r )


    def test_DrawR_r0_equal_sigma( self ):
        D = 1e-12
        kf = 1e-8
        a = 1e-7
        r0 = -a

        t = 1e-3
        
        gf = mod.FirstPassageGreensFunction1DRad( D, kf  )
        gf.seta( a )

        r = gf.drawR( 0.5, r0, t )
        self.failIf( r < -a or r > a )

    def test_DrawR_squeezed( self ):

        D = 1e-12
        kf = 1e-8
        a = 0.01e-8
        
        gf = mod.FirstPassageGreensFunction1DRad( D, kf  )
        gf.seta( a )

        t = 1e-6
        r0 = 0
        r = gf.drawR( 0.5, r0, t )
        self.failIf( r < -a or r > a )

        # near s
        r0 = -a + 0.0001e-8
        r = gf.drawR( 0.5, r0, t )
        self.failIf( r < -a or r > a )

        # near a
        r0 = a -0.0001e-8
        r = gf.drawR( 0.5, r0, t )
        self.failIf( r < -a or r > a )

        
if __name__ == "__main__":
    unittest.main()
