//#define NDEBUG
//#define BOOST_DISABLE_ASSERTS

#include <iostream>
#include <stdexcept>
#include <vector>
#include <sstream>

#include <boost/bind.hpp>
#include <boost/tuple/tuple.hpp>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_lambert.h>
#include <gsl/gsl_integration.h>

#include "factorial.hpp"
#include "funcSum.hpp"
#include "findRoot.hpp"
#include "freeFunctions.hpp"

#include "FirstPassagePairGreensFunction2D3.hpp"

const Real FirstPassagePairGreensFunction2D::TOLERANCE;
const Real FirstPassagePairGreensFunction2D::MIN_T_FACTOR;
const unsigned int FirstPassagePairGreensFunction2D::MAX_ORDER;
const unsigned int FirstPassagePairGreensFunction2D::MAX_ALPHA_SEQ;


// This is the constructor
FirstPassagePairGreensFunction2D::
FirstPassagePairGreensFunction2D( const Real D, 
				const Real kf, 
				const Real Sigma )
    :
    PairGreensFunction( D, kf, Sigma ),
    h( kf / D ),
    a( INFINITY )
//    alpha0_threshold( 0.0 )
{
    ; // do nothing
}

FirstPassagePairGreensFunction2D::~FirstPassagePairGreensFunction2D()
{
    ; // do nothing
}

void FirstPassagePairGreensFunction2D::seta( const Real a )
// This sets the parameter a (not seta or smt)
// a is the outer shell 
{
    const Real sigma( this->getSigma() );

    THROW_UNLESS( std::invalid_argument, a >= sigma );

    if( this->a != a )
    {
        this->a = a;
        clearAlphaTable();
    }
}

//
// Alpha-related methods
//

void FirstPassagePairGreensFunction2D::clearAlphaTable() const
{
    std::for_each( this->alphaTable.begin(), this->alphaTable.end(),
		   boost::mem_fn( &RealVector::clear ) );
    this->alphaOffsetTable[0] = 0;
    std::fill( this->alphaOffsetTable.begin()+1, this->alphaOffsetTable.end(),
	       -1 );

}

// The method evaluates the equation for finding the alphas for given alpha. This
// is needed to find the alpha's at which the expression is zero -> alpha is the root.
// This one is NOT USED
const Real 
FirstPassagePairGreensFunction2D::f_alpha0( const Real alpha ) const
{
	const Real a( this->geta() );
	const Real sigma( getSigma() );
	const Real h( this->geth() );
	const Real s_An( sigma * alpha );
	const Real a_An( a * alpha );

	const double J0_s_An (gsl_sf_bessel_J0(s_An));
	const double J1_s_An (gsl_sf_bessel_J1(s_An));
	const double J0_a_An (gsl_sf_bessel_J0(a_An));

	const double Y0_s_An (gsl_sf_bessel_Y0(s_An));
	const double Y1_s_An (gsl_sf_bessel_Y1(s_An));
	const double Y0_a_An (gsl_sf_bessel_Y0(a_An));

	const double rho1 ( ( (h * J0_s_An) + (alpha * J1_s_An) ) * Y0_a_An );
	const double rho2 ( ( (h * Y0_s_An) + (alpha * Y1_s_An) ) * J0_a_An );
	return rho1 - rho2;
}


// completely unclear why this is a separate function
const Real 
FirstPassagePairGreensFunction2D::f_alpha0_aux_F( const Real alpha,
						const f_alpha0_aux_params* const params )
{
    const FirstPassagePairGreensFunction2D* const gf( params->gf ); 	// get the gf from the params?
    const Real value( params->value );					// get value from the parameters

    return gf->f_alpha0( alpha );
}

// finding the ith root for the order of Bessel functions zero. Only needed for the survival probability
// and the flux
const Real 
FirstPassagePairGreensFunction2D::alpha0_i( const Integer i ) const
{
    THROW_UNLESS( std::out_of_range, i >= 0 );

    const Real a( this->geta() );		// get the a (the outer boundary
    const Real sigma( this->getSigma() );	// and sigma


    const Real target( i * M_PI + M_PI_2 );	// i*PI + PI/2
    f_alpha0_aux_params params = { this, target };	// struct with the gf object and the target


    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &FirstPassagePairGreensFunction2D::f_alpha0_aux_F ),
	    &params 
	};


    // We know the range of the solution from - Pi/2 <= atan <= Pi/2.
    const Real interval( M_PI / ( a - sigma ) );
    Real low( i * interval + std::numeric_limits<Real>::epsilon() );
    Real high( (i+1) * interval );

    //assert( GSL_FN_EVAL( &F, low ) * GSL_FN_EVAL( &F, high ) < 0.0 );

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, low, high );

    const unsigned int maxIter( 100 );

    unsigned int j( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );

        low = gsl_root_fsolver_x_lower( solver );
        high = gsl_root_fsolver_x_upper( solver );
	const int status( gsl_root_test_interval( low, high, 0.0, 1e-15 ) );

	if( status == GSL_CONTINUE )
	{
	    if( j >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "alpha0_i: failed to converge." 
			  << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}

	++j;
    }

    const Real alpha( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );
  
//    printf("alphavalue %g, alpha %g\n",GSL_FN_EVAL( &F, alpha ), alpha);

    return alpha;
}

// calculates the value of the expression for which the roots have to be found for given order n
// The expression is evaluated for given alpha, to check if it zero or to get the value of the expression
// n is the summation index used in the solution
const Real FirstPassagePairGreensFunction2D::f_alpha( const Real alpha,
						    const Integer n ) const
{
	const Real a( this->geta() );
	const Real sigma( getSigma() );
	const Real h( this->geth() );
	const Real s_An( sigma * alpha );
	const Real a_An( a * alpha );
	const Real realn( static_cast<Real>( n ) );
	const Real hs_n( h*sigma - realn );

	const double Jn_s_An  (gsl_sf_bessel_Jn(n, s_An));
	const double Jn1_s_An (gsl_sf_bessel_Jn(n+1, s_An));
	const double Jn_a_An  (gsl_sf_bessel_Jn(n, a_An));

	const double Yn_s_An  (gsl_sf_bessel_Yn(n, s_An));
	const double Yn1_s_An (gsl_sf_bessel_Yn(n+1, s_An));
	const double Yn_a_An  (gsl_sf_bessel_Yn(n, a_An));

	const double rho1 ( ( (hs_n * Jn_s_An) + (s_An * Jn1_s_An) ) * Yn_a_An );
	const double rho2 ( ( (hs_n * Yn_s_An) + (s_An * Yn1_s_An) ) * Jn_a_An );
//	cout << "hallo\n";
	return (rho1 - rho2); // /(Yn_a_An*Jn_a_An);
}


// Finding the roots alpha using an approximation for the Bessel functions
const Real 
FirstPassagePairGreensFunction2D::f_alpha_aux_F( const Real alpha,
	       					const f_alpha_aux_params* const params )
{
    const FirstPassagePairGreensFunction2D* const gf( params->gf ); 
    const Integer n( params->n );
    const Real value( params->value );


//    return gf->f_alpha_aux( alpha, n ) - value;
    return gf->f_alpha( alpha, n );	// f_alpha is the function of which we have to find the roots
}


// This calculates a root based on the previous one (offset) with Bessel functions of order n
// The roots are calculated based on the previous root which is given as an argument (offset)
const Real 
FirstPassagePairGreensFunction2D::alpha_i( const Real offset, const Integer n, 
					 gsl_root_fsolver* const solver ) const
{
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );
    const Real interval( M_PI / ( a - sigma ) );	// interval is the estimated distance to the next
							// root

    const Real target( offset + interval );	// the target root is the previous one plus the interval

    Real low ( target - (0.5 * interval ) );	// the range is calculated very easily
    Real high( target + (0.5 * interval ) );	// We are sure that there is only one root in the interval
	if (low <= 0)
	{	low = 1e-5;			// NEW to avoid the zero and negative values for the interval
	}

    f_alpha_aux_params params = { this, n, 0 };	// n is the summation index (the order of the Bessel
						// functions used
    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>
	    ( &FirstPassagePairGreensFunction2D::f_alpha_aux_F ),
	    &params 
	};

    gsl_root_fsolver_set( solver, &F, low, high );

    const unsigned int maxIter( 100 );
    unsigned int k( 0 );
    int status;

    do
    {
	gsl_root_fsolver_iterate( solver );
	
	low = gsl_root_fsolver_x_lower( solver );
	high = gsl_root_fsolver_x_upper( solver );
	status = gsl_root_test_interval( low, high, 1e-6, 1e-15 );
	
	if( status == GSL_CONTINUE && k >= maxIter )	// iteration should continue but there are too many
	{						// steps
		gsl_root_fsolver_free( solver );
		std::cerr << "alpha_i: failed to converge." 
			  << std::endl;
		throw std::exception();
	}
	++k;
    }
    while (status == GSL_CONTINUE);

    const Real alpha( gsl_root_fsolver_root( solver ) );	// get the root
    return alpha;
}


// Calculates the first root (called the offset) for summation index n (so using Bessel functions of order n)
const Real
FirstPassagePairGreensFunction2D::alphaOffset( const unsigned int n ) const
{
    if( this->alphaOffsetTable[n] > 0 )		// if the offset has already been calculated
    {
	return this->alphaOffsetTable[n];	// don't do all the work again
    }
    assert( this->alphaOffsetTable.size() >= n );	// assume the table is large enough to
							// accomodate the requested index

    const Real sigma( this->getSigma() );
    const Real a( this->geta() );
    const Real interval( M_PI_2 / ( a - sigma ));	// the x-axis of the function scales with a-sigma

	Real offset (0);	// by default the offset is zero
	if (n != 0)
	{	offset = ( this->alphaOffsetTable[n-1] );
	}						// The offsets only get bigger, the roots shift to
							// the right with the order of the Bessel functions n

    // Here we find the interval where the first positive root is in.
    // We keep shifting the search window with the size of 'interval'
    // 'interval' is half of Pi/(a - sigma) (the expected distance between the 'offsets'

    Real low ( offset );		// the searching interval for the new offset is between the previous
    Real high( offset + interval );	// one a interval further
	if (low <= 0)			// NEW to make sure things don't get negative
	{	low = 1e-5;		// we assume that the high value is never negative
	}

    Real f_low ( f_alpha(low,n) );		// get the values of the expression for the roots for these
    Real f_high( f_alpha(high,n) );		// values of alpha (at the ends of the search range)

    while( f_low * f_high > 0 ) // if the expression changes sign then there is a root in the range
    {				// and we should exit the loop
	low =  high;				// shift everything so we check out the next interval
	f_low = f_high;				// just pass the value on

	high += interval;
	f_high = f_alpha( high, n );
    }

   // We now determined roughly the interval where the first root should be.
   // Next we find the exact root in the interval.

   f_alpha_aux_params params = { this, n, 0};
   gsl_function F = 	{reinterpret_cast<typeof(F.function)>
			( &FirstPassagePairGreensFunction2D::f_alpha_aux_F ),
			&params
			};

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );  // initialize the solver
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    offset = findRoot( F, solver, low, high, 0.0,          // find the intersection between the random
                            TOLERANCE, "alphaOffset" );    // number and the cumm probability
    gsl_root_fsolver_free( solver );

    this->alphaOffsetTable[n] = offset;	// The offset found is now the first root
    return offset;
}


// calculates the constant part of the i-th term for the survival probability
const Real 
FirstPassagePairGreensFunction2D::p_survival_i( const Real alpha,
					      const Real r0 ) const
{
	const Real a( geta() );		// get the needed parameters
	const Real sigma( getSigma() );
        const Real h (this->geth());	// is this the correct one or do I need to use kf?

	const Real s_An (sigma*alpha);
	const Real a_An (a*alpha);
	const Real r0An (r0*alpha);

        const Real J0_aAn  (gsl_sf_bessel_J0(s_An));	// calculate all the required Bessel functions
        const Real J1_aAn  (gsl_sf_bessel_J1(s_An));
        const Real J0_bAn  (gsl_sf_bessel_J0(a_An));
	const Real J1_bAn  (gsl_sf_bessel_J1(a_An));

        const Real J0_r0An (gsl_sf_bessel_J0(r0An));
        const Real Y0_bAn  (gsl_sf_bessel_Y0(a_An));
        const Real Y0_r0An (gsl_sf_bessel_Y0(r0An));

	const Real Y1_bAn  (gsl_sf_bessel_Y1(a_An));
	const Real Y1_aAn  (gsl_sf_bessel_Y1(s_An));

	// calculate An,0
        const Real alpha_sq (alpha*alpha);

        const Real rho (h*J0_aAn + alpha*J1_aAn);///J0_bAn);
        const Real rho_sq (rho*rho);

        const Real B_n_0 (J0_r0An*Y0_bAn - Y0_r0An*J0_bAn);

        const Real A_i_0 ((alpha_sq * rho_sq * B_n_0)/( rho_sq - J0_bAn*J0_bAn*(h*h + alpha_sq)));

	// calculate the integral over Bn,0
	const Real B_n_0_int_tmp (Y0_bAn*( a*J1_bAn - sigma*J1_aAn ) - J0_bAn*( a*Y1_bAn - sigma*Y1_aAn ));
	const Real B_n_0_int (B_n_0_int_tmp/alpha);

	// return the total result
	const Real result (A_i_0 * B_n_0_int);
	return result;
}

// calculates the An,0 terms for determination of the flux through the outer interface
const Real 
FirstPassagePairGreensFunction2D::leavea_i( const Real alpha,
					  const Real r0 ) const
{
        const Real a( geta() );         // get the needed parameters
        const Real sigma( getSigma() );
        const Real h (this->geth());    // is this the correct one or do I need to use kf?

        const Real s_An (sigma*alpha);
        const Real a_An (a*alpha);
        const Real r0An (r0*alpha);

        const Real J0_aAn  (gsl_sf_bessel_J0(s_An));    // calculate all the required Bessel functions
        const Real J1_aAn  (gsl_sf_bessel_J1(s_An));
        const Real J0_bAn  (gsl_sf_bessel_J0(a_An));

        const Real J0_r0An (gsl_sf_bessel_J0(r0An));
        const Real Y0_bAn  (gsl_sf_bessel_Y0(a_An));
        const Real Y0_r0An (gsl_sf_bessel_Y0(r0An));

        // calculate An,0
        const Real alpha_sq (alpha*alpha);

        const Real rho (h*J0_aAn + alpha*J1_aAn);///J0_bAn);
        const Real rho_sq (rho*rho);

        const Real B_n_0 (J0_r0An*Y0_bAn - Y0_r0An*J0_bAn);

        const Real A_i_0 ((alpha_sq * rho_sq * B_n_0)/( rho_sq - J0_bAn*J0_bAn*(h*h + alpha_sq)));

	return A_i_0;
}

// Calculates the n-th term of the summation for calculating the flux through the inner interface (reaction)
const Real 
FirstPassagePairGreensFunction2D::leaves_i( const Real alpha,
					  const Real r0 ) const
{
        const Real a( geta() );         // get the needed parameters
        const Real sigma( getSigma() );
        const Real h (this->geth());    // is this the correct one or do I need to use kf?

        const Real s_An (sigma*alpha);
        const Real a_An (a*alpha);
        const Real r0An (r0*alpha);

        const Real J0_aAn  (gsl_sf_bessel_J0(s_An));    // calculate all the required Bessel functions
        const Real J1_aAn  (gsl_sf_bessel_J1(s_An));
        const Real J0_bAn  (gsl_sf_bessel_J0(a_An));
        const Real J1_bAn  (gsl_sf_bessel_J1(a_An));

        const Real J0_r0An (gsl_sf_bessel_J0(r0An));
        const Real Y0_bAn  (gsl_sf_bessel_Y0(a_An));
        const Real Y0_r0An (gsl_sf_bessel_Y0(r0An));

        const Real Y1_bAn  (gsl_sf_bessel_Y1(a_An));
        const Real Y1_aAn  (gsl_sf_bessel_Y1(s_An));

        // calculate An,0
        const Real alpha_sq (alpha*alpha);

        const Real rho (h*J0_aAn + alpha*J1_aAn);///J0_bAn);
        const Real rho_sq (rho*rho);

        const Real B_n_0 (J0_r0An*Y0_bAn - Y0_r0An*J0_bAn);

        const Real A_i_0 ((alpha_sq * rho_sq * B_n_0)/( rho_sq - J0_bAn*J0_bAn*(h*h + alpha_sq)));

	// calculate Bn,0(sigma, alpha)
	const Real B_n_0_sigma (Y0_bAn * J1_aAn - J0_bAn * Y1_aAn);

	// calculate the total result
	const Real result (alpha * A_i_0 * B_n_0_sigma);
        return result;
}


// calculates a table with all the constant factors for the survival probability
void 
FirstPassagePairGreensFunction2D::createPsurvTable( RealVector& table,
							const Real r0 ) const
{
    const RealVector& alphaTable_0( this->getAlphaTable( 0 ) );	// get the roots for the survival probability

    table.clear();				// empty the table
    table.reserve( alphaTable_0.size() );	// and get the nescessary memory

    std::transform( alphaTable_0.begin(), alphaTable_0.end(),
		    std::back_inserter( table ),
		    boost::bind( &FirstPassagePairGreensFunction2D::p_survival_i,
				 this, _1, r0 ) );	// This gets all the roots from 'begin' to 'end'
							// passes them as an argument to p_survival_i and
							// the result is passed to back_inserter
}

// Creates the tables with various Bessel functions used in drawR, the table is used to speed things up
void
FirstPassagePairGreensFunction2D::createY0J0Tables( RealVector& Y0_Table,
							RealVector& J0_Table,
							RealVector& Y0J1J0Y1_Table,
							const Real r0,
							const Real t ) const
{
	const RealVector& alphaTable_0( this->getAlphaTable( 0 ) );
							// get the roots for the survival probability
	Y0_Table.clear();				// empty the table
	J0_Table.clear();
	Y0J1J0Y1_Table.clear();

	Y0_Table.reserve( alphaTable_0.size() );	// and get the nescessary memory
	J0_Table.reserve( alphaTable_0.size() );
	Y0J1J0Y1_Table.reserve( alphaTable_0.size() );

	boost::tuple<Real,Real,Real> result;
	
	for (int count = 0; count < alphaTable_0.size(); count++)
	{	result = Y0J0J1_constants(alphaTable_0[count], t, r0);
		Y0_Table.push_back (result.get<0>());
		J0_Table.push_back (result.get<1>());
		Y0J1J0Y1_Table.push_back (result.get<2>());
	}
}


// Creates the values for in the tables Y0, J0 and Y0J1J0Y1
const boost::tuple<Real,Real,Real>
FirstPassagePairGreensFunction2D::Y0J0J1_constants ( const Real alpha,
							const Real t,
							const Real r0) const
{	const Real D(this->getD());
	const Real h(this->geth());
	const Real sigma(this->getSigma());
	const Real a(this->geta());

	const Real s_An (sigma*alpha);
	const Real a_An (a*alpha);
	const Real r0An (r0*alpha);

	const Real J0_aAn  (gsl_sf_bessel_J0(s_An));    // calculate all the required Bessel functions
	const Real J1_aAn  (gsl_sf_bessel_J1(s_An));
	const Real J0_bAn  (gsl_sf_bessel_J0(a_An));

	const Real J0_r0An (gsl_sf_bessel_J0(r0An));
	const Real Y0_bAn  (gsl_sf_bessel_Y0(a_An));
	const Real Y0_r0An (gsl_sf_bessel_Y0(r0An));

	const Real Y1_aAn  (gsl_sf_bessel_Y1(s_An));


        // calculate An,0
	const Real alpha_sq (alpha*alpha);
	const Real rho (h*J0_aAn + alpha*J1_aAn);///J0_bAn);
	const Real rho_sq (rho*rho);
	const Real B_n_0 (J0_r0An*Y0_bAn - Y0_r0An*J0_bAn);
	const Real A_i_0 ((alpha_sq * rho_sq * B_n_0)/( rho_sq - J0_bAn*J0_bAn*(h*h + alpha_sq)));

	// calculate the exponent with the time
	const Real expT( std::exp(-D*alpha_sq*t));
	// and the product
	const Real Ai0_expT (A_i_0 * expT / alpha);

	// calculate the large constant term in the intergral of Bn,0
	const Real Y0J1_J0Y1 (Y0_bAn*sigma*J1_aAn - J0_bAn*sigma*Y1_aAn);

	return boost::make_tuple (Ai0_expT*Y0_bAn, Ai0_expT*J0_bAn, Ai0_expT*Y0J1_J0Y1);
}


// calculates the ith term with exponent and time for the survival probability
const Real 
FirstPassagePairGreensFunction2D::p_survival_i_exp_table( const unsigned int i,
							const Real t,
							const Real r0,
							const RealVector& table ) const
{
    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * table[i];
}

// adds the exponential with the time to the sum. Needed for the calculation of the flux throught the outer
// interface
const Real 
FirstPassagePairGreensFunction2D::leavea_i_exp( const unsigned int i,
					      const Real t,
					      const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );
    return std::exp( - getD() * t * alpha * alpha ) * leavea_i( alpha, r0 );
}

// adds the exponential with the time to the sum. Needed for the inner interface (reaction)
const Real 
FirstPassagePairGreensFunction2D::leaves_i_exp( const unsigned int i,
					      const Real t,
					      const Real r0 ) const
{
    const Real alpha( this->getAlpha0( i ) );

    return std::exp( - getD() * t * alpha * alpha ) * leaves_i( alpha, r0 );
}


// calculates the Bossen function for a given r
const Real
FirstPassagePairGreensFunction2D::p_int_r_i_exp_table( const unsigned int i,
				                     const Real r,
                				     const RealVector& Y0_aAnTable,
               					     const RealVector& J0_aAnTable,
							const RealVector& Y0J1J0Y1Table ) const
{
	const Real alpha( this->getAlpha0( i ) );	// get the root An
	const Real r_An( r*alpha);

	const Real J1_rAn (gsl_sf_bessel_J1(r_An));
	const Real Y1_rAn (gsl_sf_bessel_Y1(r_An));

	const Real result (Y0_aAnTable[i]*r*J1_rAn - J0_aAnTable[i]*r*Y1_rAn - Y0J1J0Y1Table[i]);
	return result;
}


// This tries to guess the maximum number of n iterations it needs for calculating the survival probability
// Not really sure yet how this works

const unsigned int
FirstPassagePairGreensFunction2D::guess_maxi( const Real t ) const
{
    const unsigned int safety( 2 );

    if( t >= INFINITY )
    {
        return safety;
    }

    const Real D( getD() );
    const Real sigma( getSigma() );
    const Real a( geta() );

    const Real alpha0( getAlpha0( 0 ) );
    const Real Dt( D * t );
    const Real thr( ( exp( - Dt * alpha0 * alpha0 ) / alpha0 ) * 
                    this->TOLERANCE * 1e-1 );
    const Real thrsq( thr * thr );
    //printf("w %g %g\n", thr, 2 * Dt / thrsq );
    //printf("w0 %g\n", gsl_sf_lambert_W0( 2 * Dt / thrsq ) );

    if( thrsq <= 0.0 )
    {
        return this->MAX_ALPHA_SEQ;
    }

    const Real max_alpha( 1.0 /
                          ( sqrt( exp( gsl_sf_lambert_W0( 2 * Dt / thrsq ) ) *
                                  thrsq ) ) );
    const unsigned int 
        maxi( safety + 
              static_cast<unsigned int>( max_alpha * ( a - sigma ) / M_PI ) );

    return std::min( maxi, this->MAX_ALPHA_SEQ );
}


// Calculates the survival probability at a given time.
// This is a little wrapper for the p_survival_table so that you can easily calculate the survival probability
// at a given time
const Real 
FirstPassagePairGreensFunction2D::p_survival( const Real t,
					    const Real r0 ) const
{
    RealVector psurvTable;

    const Real p( p_survival_table( t, r0, psurvTable ) );

    return p;
}

// This actually calculates the Survival probability at time t given the particle was at r0 at time 0
// It uses the pSurvTable for efficiency (so you don't have to calculate all the constant factors all
// the time)
const Real 
FirstPassagePairGreensFunction2D::p_survival_table( const Real t,
						  const Real r0,
		  				  RealVector& psurvTable ) const
{
	Real p;
	const unsigned int maxi( guess_maxi( t ) );	// guess the maximum number of iterations required
        //printf("%d\n",maxi);
            
        if( psurvTable.size() < maxi + 1 )		// if the dimensions are good then this means
        {						// that the table is filled
        	IGNORE_RETURN getAlpha0( maxi );	// this updates the table of roots
                this->createPsurvTable( psurvTable, r0 );	// then the table is filled with data
        }

        p = funcSum_all( boost::bind( &FirstPassagePairGreensFunction2D::
                                          p_survival_i_exp_table,
                                          this,
                                          _1, t, r0, psurvTable ),
                         maxi );	// calculate the sum at time t

        return p*M_PI*M_PI_2;
}

// calculates the flux leaving through the inner interface at a given moment
const Real 
FirstPassagePairGreensFunction2D::leaves( const Real t,
					const Real r0 ) const
{
    const Real sigma(this->getSigma());

    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction2D::
                                        leaves_i_exp,
                                        this,
                                        _1, t, r0 ),
                           this->MAX_ALPHA_SEQ ) );

    return M_PI_2*sigma*p;
}

// calculates the flux leaving through the outer interface at a given moment
const Real 
FirstPassagePairGreensFunction2D::leavea( const Real t,
					const Real r0 ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction2D::
                                        leavea_i_exp,
                                        this,
                                        _1, t, r0 ),
                           this->MAX_ALPHA_SEQ ) );
    return p;
}


// calculates the sum of the sequence for drawR based upon the values in the tables and r
const Real
FirstPassagePairGreensFunction2D::p_int_r_table( const Real r,
               					const RealVector& Y0_aAnTable,
               					const RealVector& J0_aAnTable,
               					const RealVector& Y0J1J0Y1Table ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction2D::
                                        p_int_r_i_exp_table,
                                        this,
                                        _1, r, Y0_aAnTable, J0_aAnTable, Y0J1J0Y1Table ),
                           Y0_aAnTable.size() ) );
    return p*M_PI*M_PI_2;
}

// Used by drawTime
// Wrapper for p_survival_table for the interator to find the root for drawTime
const Real
FirstPassagePairGreensFunction2D::p_survival_table_F( const Real t,
                    					const p_survival_table_params* params )
{
    const FirstPassagePairGreensFunction2D* const gf( params->gf ); // the current gf (not sure why this is
								    // here)
    const Real r0( params->r0 );
    RealVector& table( params->table );		// table is empty but will be filled in p_survival_table
    const Real rnd( params->rnd );

    return rnd - gf->p_survival_table( t, r0, table );
}


// a wrapper to make p_int_r_table available to the iterator calculating the root
const Real
FirstPassagePairGreensFunction2D::p_int_r_F( const Real r,
					   const p_int_r_params* params )
{
    const FirstPassagePairGreensFunction2D* const gf( params->gf ); 
    const RealVector& Y0_aAnTable( params->Y0_aAnTable );
    const RealVector& J0_aAnTable( params->J0_aAnTable );
    const RealVector& Y0J1J0Y1Table( params->Y0J1J0Y1Table );
    const Real rnd( params->rnd );

    return gf->p_int_r_table( r, Y0_aAnTable, J0_aAnTable, Y0J1J0Y1Table) - rnd;
}


// Draws a first passage time, this could be an escape (through the outer boundary) or a reaction (through
// the inner boundary)
const Real FirstPassagePairGreensFunction2D::drawTime( const Real rnd, 
						     const Real r0 ) const
{
    const Real D( this->getD() );
    const Real sigma( this->getSigma() );
    const Real kf( this->getkf() );
    const Real a( this->geta() );

    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
    THROW_UNLESS( std::invalid_argument, sigma <= r0 && r0 <= a );

    if( r0 == a || a == sigma )		// when the particle is at the border or if the PD has no real size
    {
	return 0.0;
    }

    Real t_guess;			// where does this guess come from?
    Real dist;

/*    if( kf != 0 )			// not sure what this means
    {	*/
        dist = std::min( a - r0, r0 - sigma );	// take the shortest distance to a boundary
/*    }
    else
    {	
        dist = a - r0;
    }
*/
    t_guess = dist * dist / ( 4.0 * D );	// get some initial guess for the time, dr=sqrt(2dDt) with d
						// the dimensionality (2 in this case)
    t_guess *= .1;				// needed for the determination of the search range for the
						// zero for drawing the escape time

    const Real minT( std::min( sigma * sigma / D * this->MIN_T_FACTOR,
                               t_guess * 1e-6 ) );	// not sure yet how this works

    RealVector psurvTable;		// this is still empty as of now->not used

    p_survival_table_params params = { this, r0, psurvTable, rnd };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_survival_table_F ),
	    &params 
	};

    Real low( t_guess );	// put in a upper and lower limit (the picked time cannot be infinite!)
    Real high( t_guess );

    // adjust high and low to make sure that f( low ) and f( high ) straddle.
    const Real value( GSL_FN_EVAL( &F, t_guess ) );
    if( value < 0.0 )			// if the function is below zero at the guess the upper
    {					// boundary should be moved (passed the zero point)
        high *= 10;
        while( 1 )
        {
            const Real high_value( GSL_FN_EVAL( &F, high ) );
            
            if( high_value >= 0.0 )	// if the f(high) is positive they straddle
            {
                break;
            }

            if( fabs( high ) >= 1e10 )	// if high time is way too high forget about it
            {
                std::cerr << "Couldn't adjust high. F(" << high <<
                    ") = " << GSL_FN_EVAL( &F, high ) << "; r0 = " << r0 << 
                    ", " << dump() << std::endl;
                throw std::exception();
            }

//            printf( "drawTime: adjusting high: %g F = %g\n", 
//                    high, high_value );
            high *= 10;			// otherwise keep extending the interval
        }
    }
    else				// if the function is over zero (or at zero!) then the lower
    {					// boundary should be moved
        Real low_value_prev( value );
        low *= .1;

        while( 1 )
        {
            
            const Real low_value( GSL_FN_EVAL( &F, low ) );
            
            if( low_value <= 0.0 )
            {
                break;
            }
            
            // FIXME: (why? It isn't so bad)
            if( low <= minT ||
                ( low_value_prev - low_value ) < TOLERANCE ) // you KNOW that low_value_prev > low_value
            {
                std::cerr << "Couldn't adjust low.  Returning low (= "
                          << low << "); F(" << low <<
                    ") = " << GSL_FN_EVAL( &F, low ) << "; r0 = " << r0 << ", "
                          << dump() << std::endl;
                return low;
            }
            low_value_prev = low_value;

            //printf( "drawTime: adjusting low: %g, F = %g\n",
            //        low, low_value );
            low *= .1;
        }
    }

    // find the intersection of the cummulative survival probability and the randomly generated number
    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );	// initialize the solver
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    const Real t( findRoot( F, solver, low, high, 0.0, 		// find the intersection between the random
                            TOLERANCE, "drawTime" ) );		// number and the cumm probability

    gsl_root_fsolver_free( solver );

    return t;
}

// This determines based on the flux at a certain time, if the 'escape' was a reaction or a proper escape
const EventType
FirstPassagePairGreensFunction2D::drawEventType( const Real rnd, 
					       const Real r0,
					       const Real t ) const
{
    const Real D( this->getD() );
    const Real sigma( this->getSigma() );
    const Real kf( this->getkf() );
    const Real a( this->geta() );

    THROW_UNLESS( std::invalid_argument, 0 <= rnd && rnd < 1.0 );
    THROW_UNLESS( std::invalid_argument, sigma <= r0 && r0 < a );
    THROW_UNLESS( std::invalid_argument, t > 0.0 );

    if( kf == 0 )	// if there cannot be any flow through the radiating boundary it is always an escape
    {
        return ESCAPE;
    }
    
    // First, check if r0 is close only either to a or sigma relative
    // to Dt.  In such cases, the event type is always ESCAPE or REACTION,
    // respectively.   This avoids numerical instability in calculating
    // leavea() and/or leaves().

    // Here, use a rather large threshold for safety.
    const unsigned int H( 6 ); 				// what does this do? Why 6??
    const Real max_dist( H * sqrt( 4.0 * D * t ) );
    const Real a_dist( a - r0 );
    const Real s_dist( r0 - sigma );


    if( a_dist > max_dist )
    {
        if( s_dist < max_dist )
        {
            return REACTION;
        }
    }
    else // a_dist < max_dist
    {
        if( s_dist > max_dist )
        {
            return ESCAPE;
        }
    }

    const Real reaction( leaves( t, r0 ) );	// flux through rad boundary
    const Real escape( leavea( t, r0 ) );	// flux through abs boundary
    const Real value( reaction / ( reaction + escape ) );

    if( rnd <= value )  
    {
	return REACTION;   // leaves
    }
    else 
    {
	return ESCAPE;     // leavea
    }
}

// This draws a radius R at a given time, provided that the particle was at r0 at t=0
const Real FirstPassagePairGreensFunction2D::drawR( const Real rnd, 
						  const Real r0, 
						  const Real t ) const
{
    const Real D( this->getD() );
    const Real sigma( getSigma() );
    const Real a( this->geta() );

    THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, r0 >= sigma && r0 < a );

    if( t == 0.0 )		// if no time has passed
    {
	return r0;
    }

//	std::cout << "calculating the survival probability...\n";

    const Real psurv( p_survival( t, r0 ) );	// calculate the survival probability at this time
						// this is used as the normalization factor
						// BEWARE!!! This also produces the roots An and therefore
						// SETS THE HIGHEST INDEX -> side effect
						// VERY BAD PROGRAMMING PRACTICE!!

//	std::cout << "S = " << psurv << "\n";

    RealVector Y0_aAnTable;
    RealVector J0_aAnTable;
    RealVector Y0J1J0Y1Table;
    createY0J0Tables( Y0_aAnTable, J0_aAnTable, Y0J1J0Y1Table, r0, t);

    p_int_r_params params = { this, t, r0, Y0_aAnTable, J0_aAnTable, Y0J1J0Y1Table, rnd * psurv };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &p_int_r_F ),
	    &params 
	};

// TESTING 123
//	for (Real plaats=sigma; plaats < a; plaats += 1E-8)
//	{	std::cout << plaats << " " << GSL_FN_EVAL( &F, plaats )/psurv + rnd << "\n";
//	}

/*	std::cout << "the Y0_aAnTable has size: " << Y0_aAnTable.size() << " and values:\n";
	for (int index=0; index < Y0_aAnTable.size(); index++)
	{	std::cout << index << " " << Y0_aAnTable[index] << "\n";
	}
*/
    // adjust low and high starting from r0.
    // this is necessary to avoid root finding in the long tails where
    // numerics can be unstable.

    Real low( r0 );				// start with the initial position as the first guess
    Real high( r0 );

    const Real sqrt6Dt( sqrt( 4.0 * D * t ) );
    if( GSL_FN_EVAL( &F, r0 ) < 0.0 )
    {
        // low = r0
        unsigned int H( 3 );

        while( true )
        {
            high = r0 + H * sqrt6Dt;
            if( high > a )
            {
                if( GSL_FN_EVAL( &F, a ) < 0.0 )	// something is very wrong, this should never happen
                {
                    printf( "drawR: p_int_r_table( a ) < 0.0. returning a.\n" );
                    return a;
                }

                high = a;
                break;
            }

            const Real value( GSL_FN_EVAL( &F, high ) );
            if( value > 0.0 )
            {
                break;
            }

            ++H;
        }

    }
    else
    {
        // high = r0
        unsigned int H( 3 );

        while( true )
        {
            low = r0 - H * sqrt6Dt;
            if( low < sigma )
            {
                if( GSL_FN_EVAL( &F, sigma ) > 0.0 )
                {
                    printf( "drawR: p_int_r_table( sigma ) > 0.0. "
                            "returning sigma.\n" );
                    return sigma;
                }

                low = sigma;
                break;
            }

            const Real value( GSL_FN_EVAL( &F, low ) );
            if( value < 0.0 )
            {
                break;
            }

            ++H;
        }
    }


    // root finding by iteration.

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );

    const Real r( findRoot( F, solver, low, high, 1e-15,     // find the intersection between the random
                            TOLERANCE, "drawR" ) );          // number and the cumm probability

    gsl_root_fsolver_free( solver );

    return r;
}



// The calculates constant factor m,n for the drawing of theta. These factors are summed later.
const Real FirstPassagePairGreensFunction2D::p_m_alpha( const unsigned int n,
						      const unsigned int m,
						      const Real r,
						      const Real r0, 
						      const Real t ) const
{
	const Real sigma( this->getSigma() );
	const Real h( this->geth() );
	const Real a( this->geta() );
	const Real D( this->getD() );
//std::cout << "getting root... ";
	const Real alpha( this->getAlpha( m, n ) ); // get the n-th root using the besselfunctions of order m

//std::cout << "m n alpha " << m << " " << n << " " << alpha << std::endl;

	const Real alpha_sq( alpha * alpha );
	const Real realm( static_cast<Real>( m ) );
	const Real msq( realm * realm);
	const Real ssq( sigma * sigma);

	const Real s_Anm (sigma*alpha);
	const Real a_Anm (a*alpha);
	const Real r0Anm (r0*alpha);
	const Real r_Anm (r*alpha);

	// calculate the needed bessel functions
	const Real Jm_sAnm   (gsl_sf_bessel_Jn(m, s_Anm));
	const Real Jmp1_sAnm (gsl_sf_bessel_Jn(m+1, s_Anm));	// prime
	const Real Jm_aAnm   (gsl_sf_bessel_Jn(m, a_Anm));
	const Real Ym_aAnm   (gsl_sf_bessel_Yn(m, a_Anm));

	const Real Jm_r0Anm  (gsl_sf_bessel_Jn(m, r0Anm));
	const Real Ym_r0Anm  (gsl_sf_bessel_Yn(m, r0Anm));

	const Real Jm_rAnm   (gsl_sf_bessel_Jn(m, r_Anm));
	const Real Ym_rAnm   (gsl_sf_bessel_Yn(m, r_Anm));


	// calculating An,m
	const Real h_ma (h - realm/sigma);
	const Real rho (h_ma*Jm_sAnm + alpha*Jmp1_sAnm);///Jm_aAnm);
	const Real rho_sq (rho*rho);
	// calculating Bn,m(r')
	const Real B_n_m_r0 (Jm_r0Anm * Ym_aAnm  -  Ym_r0Anm * Jm_aAnm);

	const Real A_n_m ((alpha_sq * rho_sq * B_n_m_r0)/( rho_sq - Jm_aAnm*Jm_aAnm*(h*h + alpha_sq - msq/ssq)));


	// calculating Bn,m(r*)
	const Real B_n_m_r (Jm_rAnm * Ym_aAnm  -  Ym_rAnm * Jm_aAnm);


	// calculating the result
	const Real result( A_n_m * B_n_m_r * exp(-D*alpha_sq*t) );

	return result;
}


// This calculates the m-th constant factor for the drawTheta method. 
const Real 
FirstPassagePairGreensFunction2D::p_m( const Integer m,
				     const Real r,
				     const Real r0, 
				     const Real t ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction2D::
					p_m_alpha,
					this,
					_1, m, r, r0, t ),	// The m-th factor is a summation over n
			   MAX_ALPHA_SEQ, THETA_TOLERANCE ) );
    return p;
}

// this should make the table of constants used in the iteration for finding the root for drawTheta
// The index of the array is consistent with the index of the summation
void
FirstPassagePairGreensFunction2D::makep_mTable( RealVector& p_mTable,
					      const Real r, 
					      const Real r0, 
					      const Real t ) const
{
	const Real sigma( this->getSigma() );		// get some constants
	const Real a( this->geta() );
	const Real D(this->getD());

	p_mTable.clear();

	const Real p_0( this->p_m( 0, r, r0, t ) );	// This is the p_m where m is 0, for the denominator
	p_mTable.push_back( p_0 );			// put it in the table


	const Real p_1( this->p_m( 1, r, r0, t ) / p_0 );
	p_mTable.push_back( p_1 );			// put the first result in the table


	if( p_1 == 0 )
	{
		return;			// apparantly all the terms are zero? We are finished
	}

	const Real tolerance( THETA_TOLERANCE ); 
	const Real threshold( fabs( tolerance * p_1  ) );

	Real p_m_prev_abs( fabs( p_1 ) );	// store the previous term
	unsigned int m( 2 );
	while( true )
	{
		Real p_m( this->p_m( m, r, r0, t ) / p_0 );	// get the next term
//		std::cout << "m p_m: " << m << " " << p_m << "\n";

		if( ! std::isfinite( p_m ) )			// if the calculated value is not valid->exit
		{
			std::cerr << "makep_mTable: invalid value; " <<
				p_m << "( m= " << m << ")." << std::endl;
			break;
		}

		p_mTable.push_back( p_m );			// put the result in the table
        
//		std::cerr << m << " " << p_m << " " << threshold << std::endl;

		const Real p_m_abs( fabs( p_m ) );
							// truncate when converged enough.
		if( p_m_abs < threshold &&		// if the current term is smaller than threshold
            		p_m_prev_abs < threshold &&	// AND the previous term is also smaller
	    		p_m_abs <= p_m_prev_abs )	// AND the current term is smaller than the previous
		{
			break;
		}
	

		if( m >= this->MAX_ORDER )		// If the number of terms is too large
		{
			std::cerr << "p_m didn't converge." << std::endl;
			break;
		}
	
		++m;
		p_m_prev_abs = p_m_abs;
	}
}

// This method calculates the constants for the drawTheta method when the particle is at the boundary
const Real 
FirstPassagePairGreensFunction2D::dp_m_alpha_at_a( const unsigned int n,
						 const unsigned int m,
						 const Real r0, 
						 const Real t ) const
{
        const Real sigma( this->getSigma() );
        const Real h( this->geth() );
        const Real a( this->geta() );
        const Real D( this->getD() );

        const Real alpha( this->getAlpha( m, n ) ); // get the n-th root using the besselfunctions of order m

        const Real alpha_sq( alpha * alpha );
        const Real realm( static_cast<Real>( m ) );
        const Real msq( realm * realm);
        const Real ssq( sigma * sigma);

        const Real s_Anm (sigma*alpha);
        const Real a_Anm (a*alpha);
        const Real r0Anm (r0*alpha);

        // calculate the needed bessel functions
        const Real Jm_sAnm   (gsl_sf_bessel_Jn(m, s_Anm));
        const Real Jmp1_sAnm (gsl_sf_bessel_Jn(m+1, s_Anm));    // prime
        const Real Jm_aAnm   (gsl_sf_bessel_Jn(m, a_Anm));
        const Real Ym_aAnm   (gsl_sf_bessel_Yn(m, a_Anm));

        const Real Jm_r0Anm  (gsl_sf_bessel_Jn(m, r0Anm));
        const Real Ym_r0Anm  (gsl_sf_bessel_Yn(m, r0Anm));

        // calculating An,m
        const Real h_ma (h - realm/sigma);
        const Real rho (h_ma*Jm_sAnm + alpha*Jmp1_sAnm);///Jm_aAnm);
        const Real rho_sq (rho*rho);
        // calculating Bn,m(r')
        const Real B_n_m_r0 (Jm_r0Anm * Ym_aAnm  -  Ym_r0Anm * Jm_aAnm);

        const Real A_n_m ((alpha_sq * rho_sq * B_n_m_r0)/( rho_sq - Jm_aAnm*Jm_aAnm*(h*h + alpha_sq - msq/ssq)));

        // calculating the result
        const Real result( A_n_m * exp(-D*alpha_sq*t) );
	return result;
}

// Makes the sum over n for order m for the constants for the drawtheta Method
const Real 
FirstPassagePairGreensFunction2D::dp_m_at_a( const Integer m,
					   const Real r0, 
					   const Real t ) const
{
    const Real p( funcSum( boost::bind( &FirstPassagePairGreensFunction2D::
					dp_m_alpha_at_a,
					this,
					_1, m, r0, t ),
			   MAX_ALPHA_SEQ, THETA_TOLERANCE ) );

    return p;
}

// creates a tables of constants for drawTheta when the particle is at the edge of the domain
void
FirstPassagePairGreensFunction2D::makedp_m_at_aTable( RealVector& p_mTable,
						    const Real r0, 
						    const Real t ) const
{
        const Real sigma( this->getSigma() );           // get some constants
        const Real a( this->geta() );
        const Real D(this->getD());

        p_mTable.clear();

        const Real p_0( this->dp_m_at_a( 0, r0, t ) );// This is the p_m where m is 0, for the denominator
        p_mTable.push_back( p_0 );                       // put it in the table


        const Real p_1( this->dp_m_at_a( 1, r0, t ) / p_0 );
        p_mTable.push_back( p_1 );                      // put the first result in the table


        if( p_1 == 0 )
        {
                return;                 // apparantly all the terms are zero? We are finished
        }

        const Real tolerance( THETA_TOLERANCE );
        const Real threshold( fabs( tolerance * p_1  ) );

        Real p_m_prev_abs( fabs( p_1 ) );       // store the previous term
        unsigned int m( 2 );
        while( true )
        {
                Real p_m( this->dp_m_at_a( m, r0, t ) / p_0 );     // get the next term
//              std::cout << "m p_m: " << m << " " << p_m << "\n";
                if( ! std::isfinite( p_m ) )                    // if the calculated value is not valid->exit
                {
                        std::cerr << "makedp_m_at_aTable: invalid value; " <<
                                p_m << "( m= " << m << ")." << std::endl;
                        break;
                }

                p_mTable.push_back( p_m );                      // put the result in the table

//              std::cerr << m << " " << p_m << " " << threshold << std::endl;
                const Real p_m_abs( fabs( p_m ) );
                                                                // truncate when converged enough.
                if( p_m_abs < threshold &&              // if the current term is smaller than threshold
                        p_m_prev_abs < threshold &&     // AND the previous term is also smaller
                        p_m_abs <= p_m_prev_abs )       // AND the current term is smaller than the previous
                {
                        break;
                }
                if( m >= this->MAX_ORDER )              // If the number of terms is too large
                {
                        std::cerr << "p_m didn't converge." << std::endl;
                        break;
                }
                ++m;
                p_m_prev_abs = p_m_abs;
        }
}

// This calculates the m-th term of the summation for the drawTheta calculation
const Real 
FirstPassagePairGreensFunction2D::ip_theta_n( const unsigned int m,
						const Real theta,
						const RealVector& p_nTable) const
{
	const unsigned int m_p1 (m+1);
        return p_nTable[m_p1] * sin (m_p1*theta)/m_p1;
}


// calculates the cummulative probability of finding the particle at a certain theta
// It is used by the drawTheta method
// It uses the p_nTable for it to speed things up
const Real 
FirstPassagePairGreensFunction2D::ip_theta_table( const Real theta,
						const RealVector& p_nTable ) const
{
    const unsigned int maxm( p_nTable.size()-1 );	// get the length of the sum
							// it is shifted one because the first entry should
							// be used (m=0)

    const Real p( funcSum_all( boost::bind( &FirstPassagePairGreensFunction2D::
                                             ip_theta_n,
                                             this,
                                             _1, theta, p_nTable ),
                                maxm ) );
    return p;
}

// function to iterate when drawing the theta
const Real
FirstPassagePairGreensFunction2D::ip_theta_F( const Real theta,
					    const ip_theta_params* params )
{
    const FirstPassagePairGreensFunction2D* const gf( params->gf ); 
    const RealVector& p_nTable( params->p_nTable );	// table with useful constants
    const Real value( params->value );

    return theta/(M_PI*2) + (gf->ip_theta_table( theta, p_nTable )/M_PI) - value;
}


// This method draws a theta given a certain r and time (and intial condition of course)
const Real 
FirstPassagePairGreensFunction2D::drawTheta( const Real rnd,
					   const Real r, 
					   const Real r0, 
					   const Real t ) const
{
    const Real sigma( this->getSigma() );
    const Real a( this->geta() );

    // input parameter range checks.
    THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
    THROW_UNLESS( std::invalid_argument, sigma <= r0 && r0 < a );
    THROW_UNLESS( std::invalid_argument, sigma <= r && r <= a);
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    // t == 0 means no move.
    if( t == 0.0 )
    {
	return 0.0;
    }

    if (r == sigma)	// a reaction has occured, the angle is irrelevant
    {
	return 0.0;
    }

    //const Real Dt6_r0( sqrt( 6.0 * getD() * t ) / r0 );
    //const Real high( 5 * Dt6_r0 < 1.0 ? asin( 5 * Dt6_r0 ) : M_PI );
    //const Real mean( Dt6_r0 < 1.0 ? asin( Dt6_r0 ) : M_PI );
    //printf("theta high %g %g %g\n", high,mean, t );

    const Real high( M_PI );	// just find the theta between 0-pi, the distribution is symmetric around 0




    RealVector p_mTable;	// a table with constants to make calculations much faster

//	std::cout << "tabel maken\n";
    if( r == geta() )		// If the r is at the outer boundary
    {
//	std::cout << "deeltje op de rand\n";
	makedp_m_at_aTable( p_mTable, r0, t );	// making the table if the particle is on the outer boundary
    }
    else
    {
//	std::cout << "deeltje in het domein\n";
	makep_mTable( p_mTable, r, r0, t );	// making the table of constants for the regular case
    }
//	std::cout << "tabel gemaakt\n";

//    const Real ip_theta_pi( ip_theta_table( high, r, r0, t, p_mTable ) ); // get the normalization factor




    ip_theta_params params = { this, r, r0, t, p_mTable, rnd*0.5 };	// r, r0, t are not required

/*	for (double angle=0; angle < high; angle+= 0.1)
	{	std::cout << angle << " " << ip_theta_F(angle, &params) << "\n";
	}
*/
								// they are already in p_mTable
    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &ip_theta_F ),
	    &params 
	};

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, 0.0, high );


    const unsigned int maxIter( 100 );
    unsigned int i( 0 );
    int status;
    do				// find the root
    {
	gsl_root_fsolver_iterate( solver );
	const Real low( gsl_root_fsolver_x_lower( solver ) );
	const Real high( gsl_root_fsolver_x_upper( solver ) );
	status = gsl_root_test_interval( low, high, 1e-11, THETA_TOLERANCE );

	if( status == GSL_CONTINUE && i >= maxIter)	// if root finding should continue but the maximum
	{						// of iterations is exceeded
		gsl_root_fsolver_free( solver );	// throw an exception
		std::cerr << "drawTheta: failed to converge." << std::endl;
		throw std::exception();
	}
	++i;
    }
    while (status == GSL_CONTINUE);			// continue rootfinding
  
    //printf("%d\n", i );

    const Real theta = gsl_root_fsolver_root( solver );	// the root was found
    gsl_root_fsolver_free( solver );
    return theta;
}


//
// debug
//

const std::string FirstPassagePairGreensFunction2D::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() << ", sigma = " << this->getSigma() <<
	", a = " << this->geta() <<
	", kf = " << this->getkf() <<
	", h = " << this->geth() << std::endl;
    return ss.str();
}    
