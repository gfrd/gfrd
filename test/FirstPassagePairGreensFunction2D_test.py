#!/usr/bin/env python

__author__    = 'Koichi Takahashi <shafi@e-cell.org>'
__license__   = 'GPL'
__copyright__ = 'Copyright The Molecular Sciences Institute 2006-2007'


import unittest

import _gfrd as mod

import numpy


class FirstPassagePairGreensFunction2DTestCase( unittest.TestCase ):

    def setUp( self ):
        pass

    def tearDown( self ):
        pass
    
    def test_Instantiation( self ):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 1e-7

        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        self.failIf( gf == None )
        gf.seta( a )


    def test_DrawTime( self ):
	print "DrawT0"
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 1e-7
        r0 = 5e-8
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        self.failIf( t <= 0.0 or t >= numpy.inf )

        t = gf.drawTime( 0.0, r0 )
        self.failIf( t < 0.0 or t >= numpy.inf )

        t = gf.drawTime( 1 - 1e-16, r0 )
        self.failIf( t <= 0.0 or t >= numpy.inf )


    def test_DrawTime_a_equal_sigma( self ):
	print "DrawT1"
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = sigma
        r0 = a
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        self.assertEqual( 0.0, t )

    def test_DrawTime_a_near_sigma( self ):
	print "DrawT2"
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = sigma + sigma * 1e-6
        r0 = (a + sigma) * .5
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        self.failIf( t <= 0.0 or t >= numpy.inf )

    def test_DrawTime_r0_equal_a( self ):
	print "DrawT3"
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 1e-7
        r0 = a
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        self.assertEqual( 0.0, t )

    
    def test_DrawTime_r0_equal_sigma_kf_zero( self ):
	print "DrawT4"
        D = 1e-12
        kf = 0.0 # note this
        sigma = 1e-8
        a = 1e-7
        r0 = sigma
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        self.failIf( t < 0.0 or t >= numpy.inf )
	print "DrawT4 end"
    

    def no_test_DrawTime_r0_equal_sigma_kf_large( self ):
	print "DrawT5"
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 10e-7
        r0 = sigma + 1e-12
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        self.failIf( t < 0.0 or t >= numpy.inf )
	print "DrawT5 end"


    def test_DrawEventType( self ):
	print "DrawEvent"
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 1e-7
        r0 = 5e-8

        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
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
        sigma = 1e-8
        a = 1e-6 #sigma + sigma * 0.001
        r0 = 1.1e-8 #sigma+(a-sigma)/2

        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        t = gf.drawTime( 0.001, r0 )

        eventType = gf.drawEventType( 0.5, r0, t )
        self.failIf( eventType != 0 and eventType != 1 and eventType != 2 )

        eventType = gf.drawEventType( 0.0, r0, t )
        self.assertEqual( eventType, 0 )

        eventType = gf.drawEventType( 0.9999, r0, t )
        self.assertEqual( eventType, 1 )


    def test_DrawR( self ):
	print "DrawR"
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 1e-7
        r0 = 2e-8
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        t = 1e-3

        r = gf.drawR( 0.5, r0, t )
        self.failIf( r < sigma or r > a )

        r1 = gf.drawR( 0.0, r0, t )
        r2 = gf.drawR( 0.999999999999, r0, t )

        self.failIf( r1 < sigma or r1 > a )
        self.failIf( r2 < sigma or r2 > a )

        self.assertAlmostEqual( r1, sigma )
        self.assertAlmostEqual( r2, a )


    def test_DrawR_zerot( self ):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 1e-7
        r0 = 2e-8
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        t = 0.0

        r = gf.drawR( 0.5, r0, t )
        self.assertEqual( r0, r )


    def test_DrawR_r0_equal_sigma( self ):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 1e-7
        r0 = sigma

        t = 1e-3
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        r = gf.drawR( 0.5, r0, t )
        self.failIf( r < sigma or r > a )

    def test_DrawR_squeezed( self ):

        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 1.01e-8
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        t = 1e-6
        r0 = 1.005e-8
        r = gf.drawR( 0.5, r0, t )
        self.failIf( r < sigma or r > a )

        # near s
        r = 1.0001e-8
        r0 = 1.0001e-8
        r = gf.drawR( 0.5, r0, t )
        self.failIf( r < sigma or r > a )

        # near a
        r = 1.0099e-8
        r0 = 1.0099e-8
        r = gf.drawR( 0.5, r0, t )
        self.failIf( r < sigma or r > a )


    def test_DrawTheta( self ):
	print "DrawTheta"
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 1e-7
        r0 = 5e-8
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        t = gf.drawTime( 0.5, r0 )
        eventType = gf.drawEventType( 0.5, r0, t )
        r = gf.drawR( 0.5, r0, t )

        theta = gf.drawTheta( 0.5, r, r0, t )
        self.failIf( theta < 0.0 or theta > numpy.pi )

        theta = gf.drawTheta( 0.0, r, r0, t )
        self.failIf( theta < 0.0 or theta > numpy.pi )

        theta = gf.drawTheta( 0.999999, r, r0, t )
        self.failIf( theta < 0.0 or theta > numpy.pi )


    def test_DrawTheta_zerot( self ):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 1e-7
        r = 5e-8
        r0 = 5e-8
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        t = 0.0
        theta = gf.drawTheta( 0.5, r0, r0, t )
        self.assertEqual( 0.0, theta )

    def test_DrawTheta_smallt( self ):

        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 1e-7
        r = 2e-8
        r0 = 2e-8
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        t = 1e-7  # well this is not *very* small..
        theta = gf.drawTheta( 0.5, r, r0, t )

        self.failIf( theta < 0.0 or theta > numpy.pi )


    def test_DrawTheta_squeezed( self ):

        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 1.001e-8
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        t = 1e-8
        r = 1.0001e-8
        r0 = 1.0001e-8
        theta = gf.drawTheta( 0.5, r, r0, t )
        self.failIf( theta < 0.0 or theta > numpy.pi )

        # near s
        r = 1.00001e-8
        r0 = 1.00001e-8
        theta = gf.drawTheta( 0.5, r, r0, t )
        self.failIf( theta < 0.0 or theta > numpy.pi )

        # near a
        r = 1.00099e-8
        r0 = 1.00099e-8
        theta = gf.drawTheta( 0.5, r, r0, t )
        self.failIf( theta < 0.0 or theta > numpy.pi )


    def test_DrawTheta_r0_equal_sigma( self ):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 1e-7
        r0 = sigma

        t = 1e-3
        r = r0
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        theta = gf.drawTheta( 0.5, r, r0, t )
        self.failIf( theta < 0.0 or theta > numpy.pi )

    def test_DrawTheta_r_equal_a( self ):
        D = 1e-12
        kf = 1e-8
        sigma = 1e-8
        a = 1e-7
        r0 = 9e-8

        t = 1e-3
        r = a

        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        theta = gf.drawTheta( 0.5, r, r0, t )
        self.failIf( theta < 0.0 or theta > numpy.pi )

    def test_DrawTheta_1( self ):
        r0 =  1.0206416181e-07
        t =  4.41358538629e-08
        D = 4e-11
        sigma = 1e-07
        a = 1.05134e-07
        kf = 0 # h = 0

        r =  1.03421535312e-07
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        theta = gf.drawTheta( 0.5, r, r0, t )
        self.failIf( theta < 0.0 or theta > numpy.pi )

    def test_Alpha0( self ):

        D = 1e-12
        sigma = 1e-8
        kf = 1e-8
        
        a = 1e-7
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )
        maxerror = 0.0
        
        for i in range(100):
            alpha = gf.getAlpha0( i )
            error = abs( gf.f_alpha0( alpha ) )
            #print error/alpha, gf.f_alpha0( alpha*1.1 )/alpha
            maxerror = max( error, maxerror )

        self.failIf( abs( maxerror ) > 1e-10 )

    def test_psurvival_is_pleaves_plus_pleavea( self ):
	print "p_survival"
        D = 1e-12
        sigma = 1e-8
        kf = 1e-8

        t = 1e-5
        r0 = 5e-8
        
        a = 1e-7
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        surv = 1 - gf.p_survival( t, r0 )	## the amount of probability that left the system
        pleaves = 0.5 * t * gf.leaves( t, r0 )	## linearize the flux and calculate the integral,
        pleavea = 0.5 * t * gf.leavea( t, r0 )	## so the amount of probability left through this interface
        #print 'pll', surv, pleaves, pleavea
        self.failIf( surv <= 0.0 )
        #self.failIf( pleavea <= 0.0 or pleaves <= 0.0 )	# This fails for small times!!
        self.assertAlmostEqual( surv, pleaves + pleavea )

    def test_psurvival_smallt( self ):

        D = 1e-12
        sigma = 1e-8
        kf = 1e-8

        t = 1e-4
        r0 = 2e-8

        a = 1e-7

        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )

        for i in range( 5 ):
            psurv = 1.0 - gf.p_survival( t, r0 )	## Do the same as above.
            pleaves = 0.5 * t * gf.leaves( t, r0 ) 
            pleavea = 0.5 * t * gf.leavea( t, r0 )
            self.assertNotEqual( 1.0, psurv )
            self.assertAlmostEqual( pleaves + pleavea, psurv )
            t *= .1
    
    def test_Alphan( self ):

        D = 1e-12
        sigma = 1e-8
        kf = 1e-18
        
        a = 2e-7
        
        gf = mod.FirstPassagePairGreensFunction2D( D, kf, sigma )
        gf.seta( a )
        maxerror = 0
        
        for n in range(20):		## n cannot be so large -> limitation in order Bessel functions
            for i in range(1000):
                alpha = gf.getAlpha( n, i )
                error = abs( gf.f_alpha( alpha, n ) )
                maxerror = max( error, maxerror )

	print maxerror
        self.failIf( abs( maxerror ) > 1e-8 )
     

if __name__ == "__main__":
    unittest.main()
