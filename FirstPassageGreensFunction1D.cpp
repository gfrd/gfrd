#include <sstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <vector>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_roots.h>

#include <math.h>

#include "findRoot.hpp"
#include "FirstPassageGreensFunction1D.hpp"
#include "Defs.hpp"


// Calculates the probability of finding the particle inside the domain at 
// time t
const Real FirstPassageGreensFunction1D::p_survival (const Real t) const
{
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );

	const Real L(this->getL());
	const Real D(this->getD());
	const Real r0(this->getr0());

	if ( L < 0 || fabs (L-r0) < EPSILON || r0 < EPSILON )
	{
		// The survival probability of a zero domain is zero?
		return 0.0;
	}

	Real sum = 0, term = 0, prev_term = 0;
	Real nPI;
	const Real expo(-D*t/(L*L));	// exponent -D n^2 PI^2 t / l^2
	const Real r0_L(r0/L);
	Real n=1;

	do
	{
		if (n >= MAX_TERMEN )
		{
			std::cerr << "Too many terms for p_survival. N: " << n << std::endl;
			break;
		}

		nPI = n*M_PI;
		prev_term = term;
		term = exp(nPI*nPI*expo) * sin(nPI*r0_L) * (1.0 - cos(nPI)) / nPI;
		sum += term;
		n++;
	}
	// Is 1 a good measure or will this fail at some point?
	while (fabs(term/sum) > EPSILON*1.0 ||
		fabs(prev_term/sum) > EPSILON*1.0 ||
		n < MIN_TERMEN );

	return sum*2.0;
}

// Calculates the probability density of finding the particle at location r at 
// time t.
const Real FirstPassageGreensFunction1D::prob_r (const Real r, const Real t) const
{
	const Real L(this->getL());
	const Real r0(this->getr0());
	const Real D(this->getD());

	THROW_UNLESS( std::invalid_argument, 0 <= r && r <= L );
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );

	// if there was no time or no movement
        if (t == 0 || D == 0)
        {
		// the probability density function is a delta function
		if (r == r0)
                {
			return INFINITY;
                }
                else
                {      
			return 0.0;
                }
        }
	else if ( r < EPSILON || (L-r) < EPSILON || L < 0 )
	{	
		return 0.0;
	}

	const Real expo(-D*t/(L*L));
	const Real r0_L(r0/L);
	int n=1;
	Real nPI;
	Real sum = 0, term = 0, prev_term = 0;

	do
	{
		if (n >= MAX_TERMEN )
                {
			std::cerr << "Too many terms for prob_r. N: " << n << std::endl;
                        break;
                }

		prev_term = term;

		nPI = n*M_PI;
		term = exp(nPI*nPI*expo) * sin(nPI*r0_L) * sin(nPI*r/L);
		sum += term;
		n++;	
	}
	// Is 1E3 a good measure for the probability density?!
	while (fabs(term/sum) > EPSILON*PDENS_TYPICAL ||
		fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
		n <= MIN_TERMEN);

	return sum*2.0/L;
}

// Calculates the probability density of finding the particle at location z at 
// timepoint t, given that the particle is still in the domain.
const Real FirstPassageGreensFunction1D::calcpcum (const Real r, const Real t) const
{
	return prob_r (r, t)/p_survival (t);
}

double FirstPassageGreensFunction1D::drawT_f (double t, void *p)
{	
	// casts p to type 'struct drawT_params *'
	struct drawT_params *params = (struct drawT_params *)p;
	Real sum = 0, term = 0, prev_term = 0;
	Real Xn, exponent;
	// the maximum number of terms in the params table
	int    terms = params->terms;
	// the timescale used
	Real   tscale = params->tscale;

	int n=0;
	do
	{	if ( n >= terms )
                {       std::cerr << "Too many terms needed for DrawTime. N: " << n << std::endl;
                        break;
                }
		prev_term = term;

		Xn = params->Xn[n];
		exponent = params->exponent[n];
		term = Xn * exp(exponent * t);
		sum += term;
		n++;
	}
	while (fabs(term/sum) > EPSILON*tscale ||
                fabs(prev_term/sum) > EPSILON*tscale ||
                n <= MIN_TERMEN );

	// the intersection with the random number
	return 1.0 - 2.0*sum - params->rnd;
}

// Calculates the amount of flux leaving the left boundary at time t
const Real FirstPassageGreensFunction1D::leaves(const Real t) const
{
        THROW_UNLESS( std::invalid_argument, t >= 0.0 );

        const Real L(this->getL());
        const Real D(this->getD());
        const Real r0(this->getr0());

        if ( L < 0 || r0 < EPSILON )
	{
		// The flux of a zero domain is zero? Also if the particle 
		// started on the left boundary
                return INFINITY;
        }
	else if ( t < EPSILON*this->t_scale )
        {
		// if t=0.0 the flux must be zero
		return 0.0;
        }


        Real sum = 0, term = 0, prev_term = 0;
        Real nPI;
	const Real D_L_sq(D/(L*L));
        const Real expo(-D_L_sq*t);    // exponent -D n^2 PI^2 t / l^2
        const Real r0_L(r0/L);
        Real n=1;

        do
        {
		if (n >= MAX_TERMEN )
                {
			std::cerr << "Too many terms for p_survival. N: " << n << std::endl;
                        break;
                }

                nPI = n*M_PI;
                prev_term = term;
                term = n * exp(nPI*nPI*expo) * sin(nPI*r0_L);
                sum += term;
                n++;
        }
	// Is PDENS_TYPICAL a good measure or will this fail at some point?
        while (fabs(term/sum) > EPSILON*PDENS_TYPICAL ||
                fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
                n < MIN_TERMEN );

        return D_L_sq*2.0*M_PI*sum;
}

// Calculates the amount of flux leaving the right boundary at time t
const Real FirstPassageGreensFunction1D::leavea(const Real t) const
{
        THROW_UNLESS( std::invalid_argument, t >= 0.0 );

        const Real L(this->getL());
        const Real D(this->getD());
        const Real r0(this->getr0());

        if ( L < 0 || fabs (L-r0) < EPSILON )
        {
		// The flux of a zero domain is zero? Also if the particle 
		// started on the right boundary.
		return INFINITY;
        }
	else if ( t < EPSILON*this->t_scale )
	{
		// if t=0.0 the flux must be zero
	 	return 0.0;
	}


        Real sum = 0, term = 0, prev_term = 0;
        Real nPI;
	const Real D_L_sq(D/(L*L));
        const Real expo(-D_L_sq*t);    // exponent -D n^2 PI^2 t / l^2
        const Real r0_L(r0/L);
        Real n=1;

        do
        {       if (n >= MAX_TERMEN )
                {       std::cerr << "Too many terms for leaves. N: " << n << std::endl;
                        break;
                }

                nPI = n*M_PI;
                prev_term = term;
                term = n * exp(nPI*nPI*expo) * cos(nPI) * sin(nPI*r0_L);
                sum += term;
                n++;
        }
	// Is PDENS_TYPICAL a good measure or will this fail at some point?
        while (fabs(term/sum) > EPSILON*PDENS_TYPICAL ||
                fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
                n < MIN_TERMEN );

        return -D_L_sq*2.0*M_PI*sum;
}

// This draws an eventtype of time t based on the flux through the left (z=0) 
// and right (z=L) boundary. Although not completely accurate, it returns an 
// ESCAPE for an escape through the right boundary and a REACTION for an 
// escape through the left boundary.
const EventType FirstPassageGreensFunction1D::drawEventType( const Real rnd, const Real t ) const
{
        const Real L(this->getL());
        const Real r0(this->getr0());

        THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
	// if t=0 nothing has happened->no event!!
        THROW_UNLESS( std::invalid_argument, t > 0.0 );

        if ( fabs( r0 - L ) < EPSILON*L )
        {
		// if the particle started on the right boundary
		return ESCAPE;
        }
	else if ( r0 < EPSILON )
	{
		// if the particle started on the left boundary
		return REACTION;
	}

	const Real leaves_s (this->leaves(t));
	const Real leaves_a (this->leavea(t));
	const Real flux_total (leaves_s + leaves_a);
        const Real fluxratio (leaves_s/flux_total);

        if (rnd > fluxratio )
        {       
		return ESCAPE;
        }
        else
        {       
		return REACTION;
        }
}

// Draws the first passage time from the propensity function
const Real FirstPassageGreensFunction1D::drawTime (const Real rnd) const
{
	THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );

	const Real L(this->getL());
	const Real D(this->getD());
	const Real r0(this->getr0());

	if (D == 0.0 )
	{	
		return INFINITY;
	}
	else if (L < 0 || r0 < EPSILON || r0 > (L-EPSILON))
	{
		// if the domain had size zero
		return 0.0;
	}

	const Real expo(-D/(L*L));
	const Real r0_L(r0/L);

	// the structure to store the numbers to calculate the numbers for 1-S
	struct drawT_params parameters;
	Real nPI;
	Real Xn, exponent;
	int n = 0;

	// produce the coefficients and the terms in the exponent and put them 
	// in the params structure
	do
	{
		nPI = ((Real)(n+1))*M_PI;			//
		Xn = sin(nPI*r0_L) * (1.0 - cos(nPI)) / nPI; 
		exponent = nPI*nPI*expo;

		// store the coefficients in the structure
		parameters.Xn[n] = Xn;	
		// also store the values for the exponent
		parameters.exponent[n]=exponent;
		n++;
	}
	// modify this later to include a cutoff when changes are small
	while (n<MAX_TERMEN);

	// store the random number for the probability
	parameters.rnd = rnd;
	// store the number of terms used
	parameters.terms = MAX_TERMEN;
	parameters.tscale = this->t_scale;

//debugging
/*
for (double t=0; t<0.1; t += 0.0001)
{
	std::cout << t << " " << drawT_f (t, &parameters) << std::endl;
}
*/
	gsl_function F;
	F.function = &drawT_f;
	F.params = &parameters;


	// Find a good interval to determine the first passage time in
	const Real dist( std::min(r0, L-r0));

	// construct a guess: msd = sqrt (2*d*D*t)
	const Real t_guess( dist * dist / ( 2. * D ) );
	Real value( GSL_FN_EVAL( &F, t_guess ) );
	Real low( t_guess );
	Real high( t_guess );

	if( value < 0.0 )
	{
		// scale the interval around the guess such that the function 
		// straddles if the guess was too low
		do
		{	
			// keep increasing the upper boundary until the 
			// function straddles
			high *= 10;
			value = GSL_FN_EVAL( &F, high );

			if( fabs( high ) >= t_guess * 1e6 )
			{
				std::cerr << "Couldn't adjust high. F(" << high <<
				    ") = " << value << std::endl;
				throw std::exception();
			}
		}
		while ( value <= 0.0 );
	}
	else
	{
		// if the guess was too high initialize with 2 so the test 
		// below survives the first iteration
		Real value_prev( 2 );
		do			
		{
			if( fabs( low ) <= t_guess * 1e-6 || fabs(value-value_prev) < EPSILON*this->t_scale )
			{
				std::cerr << "Couldn't adjust low. F(" << low <<
					") = " << value << 
					" t_guess: " << t_guess << " diff: " << (value - value_prev) <<
					" value: " << value << " value_prev: " << value_prev <<
					" t_scale: " << this->t_scale <<
					std::endl;
				return low;
			}

			value_prev = value;	
			// keep decreasing the lower boundary until the 
			// function straddles
			low *= .1;
			// get the accompanying value
			value = GSL_FN_EVAL( &F, low );

		}
		while ( value >= 0.0 );
	}

	// find the intersection on the y-axis between the random number and 
	// the function
	// define a new solver type brent
	const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
	// make a new solver instance
	// incl typecast?
	gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
	const Real t( findRoot( F, solver, low, high, EPSILON*t_scale, EPSILON,
		"FirstPassageGreensFunction1D::drawTime" ) );

	// return the drawn time
	return t;
}

double FirstPassageGreensFunction1D::drawR_f (double z, void *p)
{	
	struct drawR_params *params = (struct drawR_params *)p;
	double sum = 0, term = 0, prev_term = 0;
	double S_Cn_An, n_l;
	int    terms = params->terms;

	int n=0;
	do
	{
		if (n >= terms )
                {
			std::cerr << "Too many terms for DrawR. N: " << n << std::endl;
                        break;
                }
		prev_term = term;

		S_Cn_An = params->S_Cn_An[n];
		n_l     = params->n_l[n];
		term = S_Cn_An * ( 1.0 - cos(n_l*z) );

		sum += term;
		n++;
	}
	// A lengthscale of 1 if implied here
	while (fabs(term/sum) > EPSILON ||
                fabs(prev_term/sum) > EPSILON ||
                n <= MIN_TERMEN );

	// find the intersection with the random number
	return sum - params->rnd;
}

// Draws the position of the particle at a given time, assuming that the 
// particle is still in the domain
const Real FirstPassageGreensFunction1D::drawR (const Real rnd, const Real t) const
{
	THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );

	const Real L(this->getL());
	const Real D(this->getD());
	const Real r0(this->getr0());

	// if there was no movement or the domain was zero
	if (D == 0.0 || L < 0.0 || t == 0.0)
	{	
		// scale the result back to 'normal' size
		return r0*(this->l_scale);
	}
	else
	{	// if the initial condition is at the boundary, raise an error
		// The particle can only be at the boundary in the ABOVE cases
		THROW_UNLESS( std::invalid_argument, EPSILON <= r0 && r0 <= (L-EPSILON) );
	}
	// else the normal case
	// From here on the problem is well defined


	// structure to store the numbers to calculate numbers for 1-S
	struct drawR_params parameters;
	Real S_Cn_An;
	Real nPI;
	int n=0;
	const Real S = 2.0/p_survival(t);
	const Real expo (-D*t/(L*L));
	const Real r0_L(r0/L);

//std::cout << S << std::endl;

	// produce the coefficients and the terms in the exponent and put them 
	// in the params structure
	do
	{
		nPI = ((Real)(n+1))*M_PI;			//
		S_Cn_An = S * exp(nPI*nPI*expo) * sin(nPI*r0_L) / nPI;

		// also store the values for the exponent
		parameters.S_Cn_An[n]= S_Cn_An;
		parameters.n_l[n]    = nPI/(L);
		n++;
	}
	while (n<MAX_TERMEN);

	// store the random number for the probability
	parameters.rnd = rnd ;
	// store the number of terms used
	parameters.terms = MAX_TERMEN;

	// find the intersection on the y-axis between the random number and 
	// the function
	gsl_function F;
	F.function = &drawR_f;
	F.params = &parameters;

/*
for (double x=0; x<2*a; x += 2*a/100)
{
	std::cout << x << " " << drawR_f (x, &parameters)+rnd << std::endl;
}
*/

	// define a new solver type brent
	const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent );
	// make a new solver instance
	// incl typecast?
	gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );
	const Real r( findRoot( F, solver, 0.0, L, L*EPSILON, EPSILON,
		"FirstPassageGreensFunction1D::drawR" ) );

	// return the drawn time
	return r*(this->l_scale);
}

