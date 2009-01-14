#define BOOST_AUTO_TEST_MAIN

#define BOOST_TEST_MODULE SphericalBesselGenerator

#include <boost/mpl/list.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/auto_unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "SphericalBesselGenerator.hpp"


const unsigned int maxn( 51 );

const SphericalBesselGenerator& 
generator( SphericalBesselGenerator::instance() );


const Real rel_tol( 1e-5 );
const Real abs_tol( 1e-11 );



#define CHECK_ERROR( n, z, a, b, abs_tol, rel_tol )      \
{\
    const Real abs_error( fabs( a - b ) );              \
    const Real rel_error( abs_error / fabs( a ) );\
    BOOST_CHECK_MESSAGE( abs_error < abs_tol || \
                         rel_error < rel_tol, \
                        "n " << n        \
                         << " z " << z                            \
                        << " abs error " << abs_error           \
                         << " rel error " << rel_error );\
} //                                                      \



BOOST_AUTO_TEST_CASE( testJ )
{
    const UnsignedInteger resolution( 300 );
    const Real maxz( std::max( 1000., static_cast<Real>( maxn * maxn ) ) * 2 );

    for( UnsignedInteger i( 0 ); i <= resolution; ++i )
    {
        const Real z( maxz * i / resolution );
        
        for( UnsignedInteger n( 0 ); n <= maxn; ++n )
        {
            const Real tj( generator.j( n, z ) );
            const Real j( gsl_sf_bessel_jl( n, z ) );
            
            //BOOST_CHECK_CLOSE( j, tj, TOLERANCE );
            CHECK_ERROR( n, z, j, tj, abs_tol, rel_tol );

            //printf("%d %g\n",n,z);
        }
    }

}


BOOST_AUTO_TEST_CASE( testY )
{
    const UnsignedInteger resolution( 300 );
    const Real maxz( std::max( 1000., static_cast<Real>( maxn * maxn ) ) * 2 );

    for( UnsignedInteger i( 1 ); i <= resolution; ++i )
    {
        const Real z( maxz * i / resolution );
        
        for( UnsignedInteger n( 0 ); n <= maxn; ++n )
        {
            const Real ty( generator.y( n, z ) );
            const Real y( gsl_sf_bessel_yl( n, z ) );
            
            CHECK_ERROR( n, z, y, ty, abs_tol, rel_tol );

            //printf("y %d %g\n",n,z);
        }
    }

}
