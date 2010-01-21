#!/usr/bin/env python

__author__    = 'Laurens Bossen'
__copyright__ = ''


import unittest

import _gfrd as mod

import numpy


class FirstPassageGreensFunction1DTestCase( unittest.TestCase ):

    def setUp( self ):
        pass

    def tearDown( self ):
        pass

    def test_Instantiation( self ):
        D = 1e-12
        L = 2e-7

        gf = mod.FirstPassageGreensFunction1D( D )
        self.failIf( gf == None )
        gf.setL( L )


    def test_DrawTime( self ):
        D = 1e-12
        L = 2e-7
        r0 = 5e-8

        gf = mod.FirstPassageGreensFunction1D( D )
        gf.setL( L )
        gf.setr0 ( r0 )

        t = gf.drawTime( 0.5 )
        self.failIf( t <= 0.0 or t >= numpy.inf )

        t = gf.drawTime( 0.0 )
        self.failIf( t < 0.0 or t >= numpy.inf )

        t = gf.drawTime( 1e-16 )
        self.failIf( t <= 0.0 or t >= numpy.inf )

        t = gf.drawTime( 1 - 1e-16 )
        self.failIf( t <= 0.0 or t >= numpy.inf )


    def test_DrawTime_a_equal_sigma( self ):
        D = 1e-12
        L = 0
        r0 = L

        gf = mod.FirstPassageGreensFunction1D( D )
        gf.setL( L )
        gf.setr0 ( r0 )

        t = gf.drawTime( 0.5 )
        self.assertEqual( 0.0, t )


    def test_DrawTime_a_near_sigma( self ):
        D = 1e-12
        L = 2e-14
        r0 = L/2

        gf = mod.FirstPassageGreensFunction1D( D )
        gf.setL( L )
        gf.setr0 ( r0 )

        t = gf.drawTime( 0.5 )
        self.failIf( t <= 0.0 or t >= numpy.inf )


    def test_DrawTime_r0_equal_a( self ):
        D = 1e-12
        L = 2e-7
        r0 = L

        gf = mod.FirstPassageGreensFunction1D( D )
        gf.setL( L )
        gf.setr0 ( r0 )

        t = gf.drawTime( 0.5 )
        self.assertEqual( 0.0, t )


    def test_DrawTime_r0_equal_sigma( self ):
        D = 1e-12
        L = 1e-7
        r0 = 0

        gf = mod.FirstPassageGreensFunction1D( D )
        gf.setL( L )
        gf.setr0 ( r0 )

        t = gf.drawTime( 0.5 )
        self.failIf( t < 0.0 or t >= numpy.inf )


    def test_DrawEventType( self ):
        D = 1e-12
        L = 2e-7
        r0 = L/2

        gf = mod.FirstPassageGreensFunction1D( D )
        gf.setL( L )
        gf.setr0 ( r0 )

        t = gf.drawTime( 0.5 )
        eventType = gf.drawEventType( 0.5, t )
        self.failIf( eventType != 0 and eventType != 1 and eventType != 2 )

        eventType = gf.drawEventType( 0.0, t )
        self.assertEqual( eventType, 0 )

        eventType = gf.drawEventType( 0.999999, t )
        self.assertEqual( eventType, 1 )


    def no_test_DrawEventType_smallt( self ):
        D = 1e-12
        L = 2e-6
        r0 = L/2

        gf = mod.FirstPassageGreensFunction1D( D )
        gf.setL( L )
        gf.setr0 ( r0 )

        t = gf.drawTime( 0.001 )

        eventType = gf.drawEventType( 0.5, t )
        self.failIf( eventType != 0 and eventType != 1 and eventType != 2 )

        eventType = gf.drawEventType( 0.0, t )
        self.assertEqual( eventType, 0 )

        eventType = gf.drawEventType( 0.9999, t )
        self.assertEqual( eventType, 1 )


    def test_DrawR( self ):
        D = 1e-12
        L = 2e-7
        r0 = L/2

        gf = mod.FirstPassageGreensFunction1D( D )
        gf.setL( L )
        gf.setr0 ( r0 )

        t = 1e-3

        r = gf.drawR( 0.5, t )
        self.failIf( r < 0 or r > L )

        r1 = gf.drawR( 0.0, t )
        r2 = gf.drawR( 0.999999999999, t )

        self.failIf( r1 != 0 )
        self.failIf( r2 < 0 or r2 > L )

        self.assertAlmostEqual( r1, 0 )
        self.assertAlmostEqual( r2, L )


    def test_DrawR_zerot( self ):
        D = 1e-12
        L = 1e-7
        r0 = L/2

        gf = mod.FirstPassageGreensFunction1D( D )
        gf.setL( L )
        gf.setr0 ( r0 )

        t = 0.0

        r = gf.drawR( 0.5, t )
        self.assertEqual( r0, r )


    def test_DrawR_r0_equal_sigma( self ):
        D = 1e-12
        L = 2e-7
        r0 = 0

        t = 0.0#1e-3

        gf = mod.FirstPassageGreensFunction1D( D )
        gf.setL( L )
        gf.setr0 ( r0 )

        # This raises an exception, which at this point we cannot catch
        # r = gf.drawR( 0.5, t )
        # self.failIf( r < 0 or r > L )

    def test_DrawR_squeezed( self ):

        D = 1e-12
        L = 0.02e-8

        gf = mod.FirstPassageGreensFunction1D( D )
        gf.setL( L )

        t = 1e-6
        r0 = 0
        gf.setr0 ( r0 )
        # r = gf.drawR( 0.5, t )
        # self.failIf( r < 0 or r > L )

        # near s
        r0 = 0.0001e-8
        gf.setr0 ( r0 )
        r = gf.drawR( 0.5, t )
        self.failIf( r < 0 or r > L )

        # near a
        r0 = L - 0.0001e-8
        gf.setr0 ( r0 )
        r = gf.drawR( 0.5, t )
        self.failIf( r < 0 or r > L )


if __name__ == "__main__":
    unittest.main()
