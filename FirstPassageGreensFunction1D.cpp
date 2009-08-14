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


// Berekent de kans dat het deeltje zich nog in het domein bevindt op tijdstip t, de survival probability
const Real FirstPassageGreensFunction1D::p_survival (const Real t) const
{
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );

	const Real L(this->getL());
	const Real D(this->getD());
	const Real r0(this->getr0());

	if ( L < 0 || fabs (L-r0) < EPSILON || r0 < EPSILON )
	{	return 0.0;	// The survival probability of a zero domain is zero?
	}

	Real sum = 0, term = 0;
	Real nPI;
	const Real expo(-D*t/(L*L));	// exponent -D n^2 PI^2 t / l^2
	const Real r0_L(r0/L);
	Real n=1;

	do
	{	nPI = n*M_PI;
		term = exp(nPI*nPI*expo) * sin(nPI*r0_L) * (1.0 - cos(nPI)) / nPI;
		sum += term;
		n++;
	}
	while (fabs(term/sum) > EPSILON*1.0 || n < MIN_TERMEN );
		// Is 1 a good measure or will this fail at some point?

	return sum*2.0;
}

// Berekent de kans om het deeltje op plaats x of y te vinden op tijdstip t
const Real FirstPassageGreensFunction1D::prob_r (const Real r, const Real t) const
{
	const Real L(this->getL());
	const Real r0(this->getr0());
	const Real D(this->getD());

	THROW_UNLESS( std::invalid_argument, 0 <= r && r <= L );
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );

        if (t == 0 || D == 0)   // if there was no time or no movement
        {       if (r == r0)    // the probability density function is a delta function
                {       return INFINITY;
                }
                else
                {       return 0.0;
                }
        }
	else if ( r < EPSILON || (L-r) < EPSILON || L < 0 )
	{	return 0.0;
	}

	const Real expo(-D*t/(L*L));
	const Real r0_L(r0/L);
	int n=1;
	Real nPI;
	Real sum = 0, term = 0;

	do
	{	nPI = n*M_PI;
		term = exp(nPI*nPI*expo) * sin(nPI*r0_L) * sin(nPI*r/L);
		sum += term;
		n++;	
	}
	while (fabs(term/sum) > EPSILON*PDENS_TYPICAL);	// Is 1E3 a good measure for the probability density?!

	return sum*2.0/L;
}

// Berekent de kans om het deeltje op plaats z te vinden op tijdstip t, gegeven dat het deeltje
// zich nog in het domein bevindt.
const Real FirstPassageGreensFunction1D::calcpcum (const Real r, const Real t) const
{	return prob_r (r, t)/p_survival (t);
}

double FirstPassageGreensFunction1D::drawT_f (double t, void *p)
{	struct drawT_params *params = (struct drawT_params *)p;// casts p naar type 'struct drawT_params *'
	Real sum = 0;
	Real Xn, exponent;
	int    terms = params->terms;

	for (int n=0; n<terms; n++)	// number of terms used
	{	Xn = params->Xn[n];
		exponent = params->exponent[n];
		sum += Xn * exp(exponent * t);
	}
	return 1.0 - 2.0*sum - params->rnd;		// het snijpunt vinden met het random getal
}

// Trekt een tijd uit de propensity function, een first passage time.
const Real FirstPassageGreensFunction1D::drawTime (const Real rnd) const
{
	THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );

	const Real L(this->getL());
	const Real D(this->getD());
	const Real r0(this->getr0());

	if (D == 0.0 )
	{	return INFINITY;
	}
	else if (L < 0 || r0 < EPSILON || r0 > (L-EPSILON))		// if the domain had size zero
	{	return 0.0;
	}

	const Real expo(-D/(L*L));
	const Real r0_L(r0/L);

	struct drawT_params parameters;	// the structure to store the numbers to calculate the numbers for 1-S
	Real nPI;
	Real Xn, exponent;
	int n = 0;

	// produce the coefficients and the terms in the exponent and put them in the params structure
	do
	{
		nPI = ((Real)(n+1))*M_PI;		// 
		Xn = sin(nPI*r0_L) * (1.0 - cos(nPI)) / nPI; 
		exponent = nPI*nPI*expo;

		parameters.Xn[n] = Xn;			// store the coefficients in the structure
		parameters.exponent[n]=exponent;	// also store the values for the exponent
		n++;
	}
	while (n<MAX_TERMEN);	// modify this later to include a cutoff when changes are small

	parameters.rnd = rnd;			// store the random number for the probability
	parameters.terms = MAX_TERMEN;		// store the number of terms used

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

	const Real t_guess( dist * dist / ( 2. * D ) );   // construct a guess: msd = sqrt (2*d*D*t)
	Real value( GSL_FN_EVAL( &F, t_guess ) );
	Real low( t_guess );
	Real high( t_guess );

	// scale the interval around the guess such that the function straddles
	if( value < 0.0 )		// if the guess was too low
	{
		do
		{	high *= 10;	// keep increasing the upper boundary until the function straddles
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
	else				// if the guess was too high
	{
		Real value_prev( 2 );	// initialize with 2 so the test below survives the first iteration
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
			low *= .1;	// keep decreasing the lower boundary until the function straddles
			value = GSL_FN_EVAL( &F, low );	// get the accompanying value

		}
		while ( value >= 0.0 );
	}

	// find the intersection on the y-axis between the random number and the function
	const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent ); // define a new solver type brent
	gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );  // make a new solver instance
									   // incl typecast?
	const Real t( findRoot( F, solver, low, high, EPSILON*t_scale, EPSILON,
		"FirstPassageGreensFunction1D::drawTime" ) );

	return t;				// return the drawn time
}

double FirstPassageGreensFunction1D::drawR_f (double z, void *p)
{	struct drawR_params *params = (struct drawR_params *)p;
	double sum = 0;
	double S_Cn_An, n_l;
	int    terms = params->terms;

	for (int n=0; n<terms; n++)	// number of terms used
	{	S_Cn_An = params->S_Cn_An[n];
		n_l     = params->n_l[n];
		sum = sum + S_Cn_An * ( 1.0 - cos(n_l*z) );
	}
	return sum - params->rnd;		// het snijpunt vinden met het random getal
}

// Berekent een positie gegeven dat het deeltje zich nog in het domein bevindt en er twee absorbing
// boundary conditions gelden
const Real FirstPassageGreensFunction1D::drawR (const Real rnd, const Real t) const
{
	THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );

	const Real L(this->getL());
	const Real D(this->getD());
	const Real r0(this->getr0());


	if (D == 0.0 || L < 0.0 || t == 0.0)	// if there was no movement or the domain was zero
	{	return r0*(this->l_scale);	// scale the result back to 'normal' size
	}
	else
	{	// if the initial condition is at the boundary, raise an error
		// The particle can only be at the boundary in the ABOVE cases
		THROW_UNLESS( std::invalid_argument, EPSILON <= r0 && r0 <= (L-EPSILON) );
	}
	// else the normal case
	// From here on the problem is well defined


	struct drawR_params parameters;	// structure to store the numbers to calculate numbers for 1-S
	Real S_Cn_An;
	Real nPI;
	int n=0;
	const Real S = 2.0/p_survival(t);
	const Real expo (-D*t/(L*L));
	const Real r0_L(r0/L);

//std::cout << S << std::endl;
	// produce the coefficients and the terms in the exponent and put them in the params structure
	do
	{
		nPI = ((Real)(n+1))*M_PI;			// 
		S_Cn_An = S * exp(nPI*nPI*expo) * sin(nPI*r0_L) / nPI;

		parameters.S_Cn_An[n]= S_Cn_An;		// also store the values for the exponent
		parameters.n_l[n]    = nPI/(L);
		n++;
	}
	while (n<MAX_TERMEN);

	parameters.rnd = rnd ;			// store the random number for the probability
	parameters.terms = MAX_TERMEN;		// store the number of terms used

	// find the intersection on the y-axis between the random number and the function
	gsl_function F;
	F.function = &drawR_f;
	F.params = &parameters;

/*
for (double x=0; x<2*a; x += 2*a/100)
{
	std::cout << x << " " << drawR_f (x, &parameters)+rnd << std::endl;
}
*/

	const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent ); // define a new solver type brent
	gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );  // make a new solver instance
									   // incl typecast?
	const Real r( findRoot( F, solver, 0.0, L, L*EPSILON, EPSILON,
		"FirstPassageGreensFunction1D::drawR" ) );

	return r*(this->l_scale);		// return the drawn time
}

/*
p = S(t) = intergral (0,l) v(z,t) uitrekenen
first passage tijd trekken uit 1-S(t)

p = v(z,t) uitrekenen (kans om deeltje op z te vinden op tijdstip t)
p = v(z,t)/S(t) uitrekenen (kans om deeltje op z te vinden wetende dat hij nog in het domein zit)
plaats trekken uit pc om deeltje daar te vinden op tijdstip t

q = q(t) = -dS(t)/dt uitrekenen (totale flux uit het systeem op tijdstip t)
q = q(t,0) uitrekenen (flux door boundary z=0)
relatieve flux uitrekenen
*/
