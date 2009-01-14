#include <sstream>
#include <iostream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>


#include "FreeGreensFunction.hpp"




const Real 
FreeGreensFunction::p_r( const Real r, 
                         const Real t ) const
{
    const Real D( getD() );
    const Real Dt( D * t );
    const Real Dt4( 4.0 * Dt );

    const Real Dt4Pi( Dt4 * M_PI );

    const Real term1( 1.0 / sqrt( gsl_pow_3( Dt4Pi ) ) );
    const Real term2( exp( - r * r / Dt4 ) );

    const Real jacobian( 4.0 * r * r * M_PI );

    return jacobian * term1 * term2;
}

const Real 
FreeGreensFunction::ip_r( const Real r, 
                          const Real t ) const
{
    const Real D( getD() );
    const Real Dt( D * t );
    const Real sqrtDt_r( 1.0 / sqrt( D * t ) );
    const Real sqrtPi_r( 1.0 / sqrt( M_PI ) );

    const Real term1( exp( - r * r / ( 4.0 * Dt ) ) * 
                      r * sqrtDt_r * sqrtPi_r );
    const Real term2( erf( r * 0.5 * sqrtDt_r ) );

    return term2 - term1;
}

const Real
FreeGreensFunction::ip_r_F( const Real r,
                            const ip_r_params* params )
{
    const FreeGreensFunction* const gf( params->gf ); 
    const Real t( params->t );
    const Real value( params->value );

    return gf->ip_r( r, t ) - value;
}



const Real 
FreeGreensFunction::drawR( const Real rnd, const Real t ) const
{
    // input parameter range checks.
    THROW_UNLESS( std::invalid_argument, rnd <= 1.0 && rnd >= 0.0 );
    THROW_UNLESS( std::invalid_argument, t >= 0.0 );

    // t == 0 or D == 0 means no move.
    if( t == 0.0 || getD() == 0.0 )
    {
	return 0.0;
    }

    ip_r_params params = { this, t, rnd };

    gsl_function F = 
	{
	    reinterpret_cast<typeof(F.function)>( &ip_r_F ),
	    &params 
	};

    Real max_r( 4.0 * sqrt( 6.0 * getD() * t ) );

    while( GSL_FN_EVAL( &F, max_r ) < 0.0 )
    {
        max_r *= 10;
    }

    const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
    gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
    gsl_root_fsolver_set( solver, &F, 0.0, max_r );

    const unsigned int maxIter( 100 );

    unsigned int i( 0 );
    while( true )
    {
	gsl_root_fsolver_iterate( solver );
	const Real low( gsl_root_fsolver_x_lower( solver ) );
	const Real high( gsl_root_fsolver_x_upper( solver ) );
	const int status( gsl_root_test_interval( low, high, 1e-15, 
						  this->TOLERANCE ) );

	if( status == GSL_CONTINUE )
	{
	    if( i >= maxIter )
	    {
		gsl_root_fsolver_free( solver );
		std::cerr << "drawR: failed to converge." << std::endl;
		throw std::exception();
	    }
	}
	else
	{
	    break;
	}

	++i;
    }
  
    //printf("%d\n", i );

    const Real r( gsl_root_fsolver_root( solver ) );
    gsl_root_fsolver_free( solver );
    
    return r;
}


const std::string FreeGreensFunction::dump() const
{
    std::ostringstream ss;
    ss << "D = " << this->getD() << std::endl;
    return ss.str();
}    
