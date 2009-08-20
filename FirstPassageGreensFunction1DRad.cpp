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

#include "FirstPassageGreensFunction1DRad.hpp"

double FirstPassageGreensFunction1DRad::tan_f (double x, void *p)
		// this is the appropriate definition of the function in gsl
{	struct tan_f_params *params = (struct tan_f_params *)p;	//  casts the void naar struct pointer
	const Real L = (params->L);
	const Real h = (params->h);
	const Real h_L (h*L);
	if ( h_L < 1 )
	{	return  1/tan(x) + (h_L)/x;		// h = k/D, h=h1/k1
	}
	else
	{	return  tan(x) + x/(h_L);		// h = k/D, h=h1/k1
	}
}

// Berekent de wortels van tan(aL)=-a*D/k
const Real FirstPassageGreensFunction1DRad::a_n(const int n) const
{
	const Real L(this->getL());			// L=length of domain=2*a
	const Real h(this->getk()/this->getD());	// h=k/D
	Real upper, lower;

	if ( h*L < 1 )
	{	lower = (n-1)*M_PI + 1E-10;	// 0.001 om ervoor te zorgen dat hij niet precies
		upper =  n   *M_PI - 1E-10;	// de overgang van -oneindig naar +oneindig meeneemt
	}
	else
	{	lower = (n-1)*M_PI + M_PI_2 + 1E-10;
		upper = n    *M_PI + M_PI_2 - 1E-10;
	}

	gsl_function F;
	struct tan_f_params params = { L, h};
     
	F.function = &FirstPassageGreensFunction1DRad::tan_f;
	F.params = &params;


	const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent ); // define a new solver type brent
	gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );  // make a new solver instance
									   // incl typecast?
	const Real a( findRoot( F, solver, lower, upper, 1.0*EPSILON, EPSILON,
		"FirstPassageGreensFunction1DRad::root_tan" ) );
	gsl_root_fsolver_free( solver );
	return a;
}

// r0 is here still in the domain from 0 - L
// The root An is also from this domain
const Real FirstPassageGreensFunction1DRad::An (const Real a_n) const
{
	const Real h(this->getk()/this->getD());
	const Real L(this->getL());
	const Real r0(this->getr0());
	const Real Anr0 = a_n*r0;

	return (a_n*cos(Anr0) + h*sin(Anr0)) / (L*(a_n*a_n + h*h) + h);
}

// r is here still in the domain from 0 - L
// The root An is also from this domain
const Real FirstPassageGreensFunction1DRad::Bn (const Real a_n) const
{	const Real h(this->getk()/this->getD());
	const Real L(this->getL());
	const Real Anr = a_n*L;
	const Real hAn = h/a_n;

	return sin(Anr) - hAn*cos(Anr) + hAn;
}

// The root An is from the domain from 0 - L
const Real FirstPassageGreensFunction1DRad::Cn (const Real a_n, const Real t) const
{	const Real D(this->getD());

	return std::exp(-D*a_n*a_n*t);
}

// Berekent de kans dat het deeltje zich nog in het domein bevindt op tijdstip t, de survival probability
// Het domein is van -r to r (r0 zit hier dus tussen!!)
const Real FirstPassageGreensFunction1DRad::p_survival (const Real t) const
{
	const Real D(this->getD());

	THROW_UNLESS( std::invalid_argument, t >= 0.0 );

	if (t == 0 || D == 0)
        {       return 1.0;      // if there was no time or no movement the particle was always in the domain
        }


	Real An;
	Real sum = 0, term = 0, term_prev = 0;
	int n = 1;

	do
	{	An = this->a_n(n);
		term_prev = term;
		term = Cn(An, t) * this->An(An) * Bn(An);
		sum += term;
		n++;
	}
	while ( fabs(term/sum) > EPSILON*1.0 ||
		fabs(term_prev/sum) > EPSILON*1.0 ||
		n <= MIN_TERMEN);
		// Is 1.0 a good measure for the scale of probability or will this fail at some point?

	return 2.0*sum;
}


// Berekent de kansdichtheid om het deeltje op plaats z te vinden op tijdstip t
const Real FirstPassageGreensFunction1DRad::prob_r (const Real r, const Real t) const
{
	const Real L(this->getL());
	const Real D(this->getD());
	const Real h(this->getk()/D);
	const Real r0(this->getr0());

	THROW_UNLESS( std::invalid_argument, t >= 0.0 );
	THROW_UNLESS( std::invalid_argument, 0 <= r && r <= L);

	if (t == 0 || D == 0)	// if there was no time or no movement
	{	if (r == r0)	// the probability density function is a delta function
		{	return INFINITY;
		}
		else
		{	return 0.0;
		}
	}

	if ( fabs (r - L) < EPSILON*L )	// if you're looking on the boundary
	{	return 0.0;
	}

	Real root_n, An_r;
	Real sum = 0, term = 0, prev_term = 0;
	int n=1;

	do
	{	if ( n >= MAX_TERMEN )
                {       std::cerr << "Too many terms needed for prob_r. N: " << n << std::endl;
                        break;
                }

		root_n = this->a_n(n);
		An_r = root_n*r;

		prev_term = term;
		term = Cn(root_n, t) * An(root_n) * (root_n*cos(An_r) + h*sin(An_r));
		sum += term;

		n++;
	}
	while (fabs(term/sum) > EPSILON*PDENS_TYPICAL || 
		fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
		n <= MIN_TERMEN );
		// PDENS_TYPICAL is now 1e3, is this any good?!

	return 2.0*sum;
}

// Berekent de kans om het deeltje op plaats z te vinden op tijdstip t, gegeven dat het deeltje
// zich nog in het domein bevindt.
const Real FirstPassageGreensFunction1DRad::calcpcum (const Real r, const Real t) const
{
	const Real r_corr(r/this->l_scale);		// BEWARE: HERE THERE IS SCALING OF R!
	return prob_r (r_corr, t)/p_survival (t);
}

// Berekent de totale kansflux die het systeem verlaat op tijdstip t
const Real FirstPassageGreensFunction1DRad::flux_tot (const Real t) const
{	Real An;
	double sum = 0, term = 0, prev_term = 0;
	const double D(this->getD());

	int n=1;

	do
	{	if ( n >= MAX_TERMEN )
                {       std::cerr << "Too many terms needed for flux_tot. N: " << n << std::endl;
                        break;
                }

		An = this->a_n(n);
		prev_term = term;
		term = An * An * Cn(An, t) * this->An(An) * Bn(An);
		n++;
		sum += term;
	}
	while (fabs(term/sum) > EPSILON*PDENS_TYPICAL ||
		fabs(prev_term/sum) > EPSILON*PDENS_TYPICAL ||
		n <= MIN_TERMEN );

	return sum*2.0*D;
}

// Berekent de kansflux die het systeem verlaat door de radiating boundary condition op x=0
// op tijdstip t
const Real FirstPassageGreensFunction1DRad::flux_rad (const Real t) const
{
	return k*prob_r(0, t);
}

// Berekent de flux door de radiating boundary condition als fractie van de totale flux
// Hiermee kunnen we de relatieve kans bepalen dat het deeltje het systeem heeft verlaten door
// de radiating boundary en niet door de absorbing boundary.
const Real FirstPassageGreensFunction1DRad::fluxRatioRadTot (const Real t) const
{	return flux_rad (t)/flux_tot (t);
}

const EventType
FirstPassageGreensFunction1DRad::drawEventType( const Real rnd, const Real t ) const
{
	const Real L(this->getL());
	const Real r0(this->getr0());

	THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
        THROW_UNLESS( std::invalid_argument, t > 0.0 );		// if t=0 nothing has happened->no event!!

	if ( k == 0 || fabs( r0 - L ) < EPSILON*L )
	{	return ESCAPE;
	}

	const Real fluxratio (this->fluxRatioRadTot(t));

	if (rnd > fluxratio )
	{	return ESCAPE;
	}
	else
	{	return REACTION;
	}
}

double FirstPassageGreensFunction1DRad::drawT_f (double t, void *p)
{	struct drawT_params *params = (struct drawT_params *)p;// casts p naar type 'struct drawT_params *'
	Real sum = 0, term = 0, prev_term = 0;
	Real Xn, exponent;
	int terms = params->terms;
	Real tscale = params->tscale;

	int n=0;
	do
	{	if ( n >= terms )
		{	std::cerr << "Too many terms needed for DrawTime. N: " << n << std::endl;
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

	return 1.0 - 2.0*sum - params->rnd;		// het snijpunt vinden met het random getal
}

// Trekt een tijd uit de propensity function, een first passage time.
const Real FirstPassageGreensFunction1DRad::drawTime (const Real rnd) const
{
	const Real L(this->getL());
	const Real k(this->getk());
	const Real D(this->getD());
	const Real r0(this->getr0());

	THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );

std::cerr << "D: " << D << " L: " << this->l_scale << " r0: " << r0 << " h: " << k/D << std::endl;

	if ( D == 0.0 || L == INFINITY )
	{	return INFINITY;
	}

	if ( rnd <= EPSILON || L < 0.0 || fabs(r0 - L) < EPSILON*L )
	{	return 0.0;
	}


	const Real h(k/D);

	struct drawT_params parameters;	// the structure to store the numbers to calculate the numbers for 1-S
	double An = 0;
	double tmp0, tmp1, tmp2, tmp3;
	double Xn, exponent;


	// produce the coefficients and the terms in the exponent and put them in the params structure
	for (int n=0; n<MAX_TERMEN; n++)
	{
		An = a_n (n+1);			// get the n-th root of tan(alfa*L)=alfa/-k
		tmp0 = An * An;			// An^2
		tmp1 = An * r0;			// An * z'
		tmp2 = An * L;			// An * L
		tmp3 = h / An;			// h / An
		Xn = (An*cos(tmp1) + h*sin(tmp1))*(sin(tmp2)-tmp3*cos(tmp2)+tmp3)/(L*(tmp0+h*h)+h); 
		exponent = -D*tmp0;

		parameters.Xn[n] = Xn;			// store the coefficients in the structure
		parameters.exponent[n]=exponent;	// also store the values for the exponent

	}
	parameters.rnd = rnd;			// store the random number for the probability
	parameters.terms = MAX_TERMEN;		// store the number of terms used
	parameters.tscale = this->t_scale;

	// Define the function for the rootfinder
	gsl_function F;
	F.function = &FirstPassageGreensFunction1DRad::drawT_f;
	F.params = &parameters;


        // Find a good interval to determine the first passage time in
        const Real dist(L-r0);			// get the distance to absorbing boundary (disregard rad BC)
        const Real t_guess( dist * dist / ( 2. * D ) );   // construct a guess: msd = sqrt (2*d*D*t)
        Real value( GSL_FN_EVAL( &F, t_guess ) );
        Real low( t_guess );
        Real high( t_guess );

//std::cout << " t_guess: " << t_guess <<
//	     " value: " << value <<
//	     std::endl;

        // scale the interval around the guess such that the function straddles
        if( value < 0.0 )               // if the guess was too low
        {
                do
                {       high *= 10;     // keep increasing the upper boundary until the function straddles
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
        else                            // if the guess was too high
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
					" rnd: " << rnd <<
					std::endl;
                                return low;
                        }
                        value_prev = value;
		        low *= .1;      // keep decreasing the lower boundary until the function straddles
                        value = GSL_FN_EVAL( &F, low );     // get the accompanying value
                }
                while ( value >= 0.0 );
        }


	// find the intersection on the y-axis between the random number and the function
	const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent ); // define a new solver type brent
	gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );  // make a new solver instance
									   // incl typecast?
	const Real t( findRoot( F, solver, low, high, t_scale*EPSILON, EPSILON,
		"FirstPassageGreensFunction1DRad::drawTime" ) );

	return t;				// return the drawn time
}

double FirstPassageGreensFunction1DRad::drawR_f (double z, void *p)
{	struct drawR_params *params = (struct drawR_params *)p;// casts p naar type 'struct drawR_params *'
	Real sum = 0, term = 0, prev_term = 0;
	Real An, S_Cn_An, b_An;
	int    terms = params->terms;

	int n = 0;
	do
	{	if ( n >= terms )
                {       std::cerr << "Too many terms needed for DrawR. N: " << n << std::endl;
                        break;
                }
		prev_term = term;

		S_Cn_An = params->S_Cn_An[n];
		b_An    = params->b_An[n];
		An      = params->An[n];
		term = S_Cn_An * (sin(An*z) -b_An * cos(An*z) + b_An);

		sum += term;
		n++;
	}
	while (fabs(term/sum) > EPSILON ||		// the lengthscale of 1.0 is implied
                fabs(prev_term/sum) > EPSILON ||
                n <= MIN_TERMEN );

	return sum - params->rnd;		// het snijpunt vinden met het random getal
}

const Real FirstPassageGreensFunction1DRad::drawR (const Real rnd, const Real t) const
{
	const Real L(this->getL());
	const Real D(this->getD());
	const Real r0(this->getr0());

        THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
        THROW_UNLESS( std::invalid_argument, t >= 0.0 );

	if (t == 0.0 || D == 0)
	{	return r0*this->l_scale;	// the trivial case
	}
	if ( L < 0.0 )			// if the domain had zero size
	{	return 0.0;
	}

	struct drawR_params parameters;	// the structure to store the numbers to calculate the numbers for 1-S
	double An = 0;
	double S_Cn_An;
	double tmp0, tmp1;
	const Real h(this->getk()/D);
	const Real S = 2.0/p_survival(t);


	// produce the coefficients and the terms in the exponent and put them in the params structure
	for (int n=0; n<MAX_TERMEN; n++)
	{
		An = a_n (n+1);			// get the n-th root of tan(alfa*L)=alfa/-k
		tmp0 = An * An;			// An^2
		tmp1 = An * r0;			// An * z'
		S_Cn_An = S * exp(-D*tmp0*t) * (An*cos(tmp1) + h*sin(tmp1)) / (L*(tmp0 + h*h) + h);

		parameters.An[n]     = An;	// store the coefficients in the structure
		parameters.S_Cn_An[n]= S_Cn_An;	// also store the values for the exponent
		parameters.b_An[n]   = h/An;

	}
	parameters.rnd = rnd;			// store the random number for the probability
	parameters.terms = MAX_TERMEN;		// store the number of terms used


	// find the intersection on the y-axis between the random number and the function
	gsl_function F;
	F.function = &FirstPassageGreensFunction1DRad::drawR_f;
	F.params = &parameters;


	const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent ); // define a new solver type brent
	gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );  // make a new solver instance
									   // incl typecast?
	const Real z( findRoot( F, solver, 0.0, L, EPSILON*L, EPSILON,
		"FirstPassageGreensFunction1DRad::drawR" ) );

	return z*this->l_scale;				// return the drawn place
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
