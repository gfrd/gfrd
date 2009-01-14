#if !defined( __FREEPAIRGREENSFUNCTION )
#define __FREEPAIRGREENSFUNCTION 

#include <cmath>

#include <gsl/gsl_integration.h>

#include "PairGreensFunction.hpp"

/**
   Pair Green's function for the case where the pair never interact.

   Therefore, drawTime() always returns +INFINITY.
   kf == sigma == 0.
*/

class FreePairGreensFunction
    :
    public PairGreensFunction
{
    
private:

    static const Real TOLERANCE = 1e-8;

    static const Real H = 7;
    
public:
    
    FreePairGreensFunction( const Real D )
        :
        PairGreensFunction( D, 0.0, 0.0 )
    {
        ; // do nothing
    }
    
    
    ~FreePairGreensFunction()
    {
        ; // do nothing
    }
    
    
    const Real drawTime( const Real rnd, const Real r0,
			 const Real maxt ) const
    {
        return INFINITY;
    }
    
    const Real drawR( const Real rnd, 
		      const Real r0, 
		      const Real t ) const;
    
    const Real drawTheta( const Real rnd,
			  const Real r, 
			  const Real r0, 
			  const Real t ) const;

    const Real p_r( const Real r, 
                    const Real r0, 
                    const Real t ) const;

    const Real ip_r( const Real r, 
                     const Real r0, 
                     const Real t ) const;
    

    const Real p_theta( const Real theta,
                        const Real r, 
                        const Real r0, 
                        const Real t ) const;

    const Real ip_theta( const Real theta,
                         const Real r, 
                         const Real r0, 
                         const Real t ) const;


    const std::string dump() const;
    

private:

    struct ip_r_params
    { 
	const FreePairGreensFunction* const gf;
	const Real r0;
	const Real t;
	const Real value;
    };


    static const Real ip_r_F( const Real r,
                              const ip_r_params* params );


    struct ip_theta_params
    { 
	const FreePairGreensFunction* const gf;
	const Real r;
	const Real r0;
	const Real t;
	const Real value;
    };


    static const Real ip_theta_F( const Real theta,
                                  const ip_theta_params* params );
    
};



#endif // __PLAINPAIRGREENSFUNCTION 
