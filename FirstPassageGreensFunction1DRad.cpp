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

#define TERMEN 500
#define THRESHOLD 1E-9

#include "FirstPassageGreensFunction1DRad.hpp"

double FirstPassageGreensFunction1DRad::tan_f (double x, void *p)
		// this is the appropriate definition of the function in gsl
{	struct tan_f_params *params = (struct tan_f_params *)p;	//  casts the void naar struct pointer
	const Real L = (params->L);
	const Real h = (params->h);
	return  1/tan(x) + (h*L)/x;		// h = k/D, h=h1/k1
}

// Berekent de wortels van tan(aL)=-a*D/k
const Real FirstPassageGreensFunction1DRad::a_n(const int n) const
{
	const Real L(this->geta()*2);			// L=length of domain=2*a
	const Real h(this->getk()/this->getD());	// h=k/D
	double lower = (n-1)*M_PI + 1E-10;	// 0.001 om ervoor te zorgen dat hij niet precies
	double upper =  n   *M_PI - 0.0001;		// de overgang van -oneindig naar +oneindig meeneemt

//	std::cout << "lower is " << lower << std::endl;
//	std::cout << "upper is " << upper << std::endl;

	gsl_function F;
	struct tan_f_params params = { L, h};
     
	F.function = &FirstPassageGreensFunction1DRad::tan_f;
	F.params = &params;
/*
for (double x=lower; x<upper; x+=((upper-lower)/100))
{
	double a = tan_f (x, &params);
	std::cout << x << " " << a << std::endl;
}
*/
	const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent ); // define a new solver type brent
	gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );  // make a new solver instance
									   // incl typecast?
	const Real a( findRoot( F, solver, lower, upper, 1e-18, 1e-12,
		"FirstPassageGreensFunction1DRad::root_tan" ) );
	gsl_root_fsolver_free( solver );
	return a/L;
}

// r0 is here still in the domain from 0 - L
// The root An is also from this domain
const Real FirstPassageGreensFunction1DRad::An (const Real a_n, const Real r0) const
{
	const Real h(this->getk()/this->getD());
	const Real L(this->geta()*2);
	const Real Anr0 = a_n*r0;

	return (a_n*cos(Anr0) + h*sin(Anr0)) / (L*(a_n*a_n + h*h) + h);
}

// r is here still in the domain from 0 - L
// The root An is also from this domain
const Real FirstPassageGreensFunction1DRad::Bn (const Real a_n, const Real r) const
{	const Real h(this->getk()/this->getD());
	const Real Anr = a_n*r;
	const Real hAn = h/a_n;

	return sin(Anr) - hAn*cos(Anr) + hAn;
}

// The root An is from the domain from 0 - L
const Real FirstPassageGreensFunction1DRad::Cn (const Real a_n, const Real t) const
{	const Real D(this->getD());

	return exp(-D*a_n*a_n*t);
}

// Berekent de kans dat het deeltje zich nog in het domein bevindt op tijdstip t, de survival probability
// Het domein is van -r to r (r0 zit hier dus tussen!!)
const Real FirstPassageGreensFunction1DRad::p_survival (const Real t, const Real r0) const
{
	const Real D(this->getD());
	const Real a(this->geta());
	THROW_UNLESS( std::invalid_argument, -a <= r0 && r0 <= a);
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );

	if (t == 0 || D == 0)
        {       return 1.0;      // if there was no time or no movement the particle was always in the domain
        }


	Real An;
	Real sum = 0, term = 0;
	int n = 1;
	const Real r0_2 (r0 + a);	// shift the coordinates with a to make it from {-a, a} to {0, L}

	do
	{	An = this->a_n(n);
		term = Cn(An, t) * this->An(An, r0_2) * Bn(An, a+a);
		sum += term;
		n++;
	}
	while (fabs(term/sum) > THRESHOLD);

	return 2.0*sum;
}


// Berekent de kansdichtheid om het deeltje op plaats z te vinden op tijdstip t
const Real FirstPassageGreensFunction1DRad::prob_r (const Real r0, const Real r, const Real t) const
{
	const Real a(geta());
	const Real D(getD());
	const Real h(this->getk()/this->getD());

	THROW_UNLESS( std::invalid_argument, t >= 0.0 );
	THROW_UNLESS( std::invalid_argument, -a <= r && r <= a);
	THROW_UNLESS( std::invalid_argument, -a <= r0 && r0 <= a);

	if (t == 0 || D == 0)	// if there was no time or no movement
	{	if (r0 == r)	// the probability density function is a delta function
		{	return INFINITY;
		}
		else
		{	return 0;
		}
	}

	Real root_n, An_r;
	Real sum = 0, term = 0;
	int n=1;
	const Real r0_2 (r0 + a);	// shift the coordinates with a to make it from {-a, a} to {0, L}
	const Real r_2  (r + a);

	do
	{	root_n = this->a_n(n);
		An_r = root_n*r_2;
		term = Cn(root_n, t) * An(root_n, r0_2) * (root_n*cos(An_r) + h*sin(An_r));
		sum += term;
		n++;
	}
	while (fabs(term/sum) > THRESHOLD);

	return 2*sum;
}

// Berekent de kans om het deeltje op plaats z te vinden op tijdstip t, gegeven dat het deeltje
// zich nog in het domein bevindt.
const Real FirstPassageGreensFunction1DRad::calcpcum (const Real r0, const Real r, const Real t) const
{	return prob_r (r0, r, t)/p_survival (t, r0);
}

// Berekent de totale kansflux die het systeem verlaat op tijdstip t
const Real FirstPassageGreensFunction1DRad::flux_tot (const Real r0, const Real t) const
{	Real An;
	double sum = 0, term = 0;
	const double D(this->getD());
	const Real a(this->geta());
	const Real L(a+a);

	int n=1;
	const Real r0_2 (r0 + a);	// shift the coordinates with a to make it from {-a, a} to {0, L}

	do
	{	An = this->a_n(n);
		term = An * An * Cn(An, t) * this->An(An, r0_2) * Bn(An, L);
		n++;
		sum += term;
	}
	while (fabs(term/sum) > THRESHOLD);

	return sum*2*D;
}

// Berekent de kansflux die het systeem verlaat door de radiating boundary condition op x=0
// op tijdstip t
const Real FirstPassageGreensFunction1DRad::flux_rad (const Real r0, const Real t) const
{	const Real D(this->getD());
	const Real k(this->getk());
	const Real a(this->geta());
	return k*prob_r(r0, -a, t);

//	return D*h*prob_r(r0, 0, t)/k;	// this is the old one, not sure if this is correct
}

// Berekent de flux door de radiating boundary condition als fractie van de totale flux
// Hiermee kunnen we de relatieve kans bepalen dat het deeltje het systeem heeft verlaten door
// de radiating boundary en niet door de absorbing boundary.
const Real FirstPassageGreensFunction1DRad::fluxRatioRadTot (const Real r0, const Real t) const
{	return flux_rad (r0, t)/flux_tot (r0, t);
}

const EventType
FirstPassageGreensFunction1DRad::drawEventType( const Real rnd, const Real r0, const Real t ) const
{
	const Real a(this->geta());

	THROW_UNLESS( std::invalid_argument, rnd < 1.0 && rnd >= 0.0 );
        THROW_UNLESS( std::invalid_argument, -a <= r0 && r0 <= a );
        THROW_UNLESS( std::invalid_argument, t > 0.0 );		// if t=0 nothing has happened->no event!!

	const Real fluxratio (this->fluxRatioRadTot(r0, t));

	if (rnd > fluxratio )
	{	return ESCAPE;
	}
	else
	{	return REACTION;
	}
}

double FirstPassageGreensFunction1DRad::drawT_f (double t, void *p)
{	struct drawT_params *params = (struct drawT_params *)p;// casts p naar type 'struct drawT_params *'
	double sum = 0;
	double Xn, exponent;
	int    terms = params->terms;

	for (int n=0; n<terms; n++)	// number of terms used
	{	Xn = params->Xn[n];
		exponent = params->exponent[n];
		sum += Xn * exp(exponent * t);
	}
//	std::cout << ".";
	return 1 - 2*sum - params->rnd;		// het snijpunt vinden met het random getal
}

// Trekt een tijd uit de propensity function, een first passage time.
const Real FirstPassageGreensFunction1DRad::drawTime (const Real rnd, const Real r0) const
{
	const Real a(this->geta());

	THROW_UNLESS( std::invalid_argument, -a <= r0 && r0 <= a );
	THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );

	const Real L(a*2);
	const Real h(this->getk()/this->getD());
	const Real r0_2(r0 + a);

	struct drawT_params parameters;	// the structure to store the numbers to calculate the numbers for 1-S
	double An = 0;
	double tmp0, tmp1, tmp2, tmp3;
	double Xn, exponent;

//std::cout << "making terms\n"; 

	// produce the coefficients and the terms in the exponent and put them in the params structure
	for (int n=0; n<TERMEN; n++)
	{
		An = a_n (n+1);		// get the n-th root of tan(alfa*L)=alfa/-k
		tmp0 = An * An;			// An^2
		tmp1 = An * r0_2;		// An * z'
		tmp2 = An * L;			// An * L
		tmp3 = h / An;			// h / An
		Xn = (An*cos(tmp1) + h*sin(tmp1))*(sin(tmp2)-tmp3*cos(tmp2)+tmp3)/(L*(tmp0+h*h)+h); 
		exponent = -D*tmp0;

		parameters.Xn[n] = Xn;			// store the coefficients in the structure
		parameters.exponent[n]=exponent;	// also store the values for the exponent

//	cout << n << " " << Xn * exp(exponent * 0.01) << endl;
	}
	parameters.rnd = rnd;			// store the random number for the probability
	parameters.terms = TERMEN;		// store the number of terms used
//std::cout << "terms made\n"; 

	// Define the function for the rootfinder
	gsl_function F;
	F.function = &FirstPassageGreensFunction1DRad::drawT_f;
	F.params = &parameters;


        // Find a good interval to determine the first passage time in
        const Real t_guess( a * a / ( 2. * D ) );   // construct a guess: msd = sqrt (2*d*D*t)
        Real value( GSL_FN_EVAL( &F, t_guess ) );
        Real low( t_guess );
        Real high( t_guess );

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
                Real value_prev( value );
                do
                {       low *= .1;      // keep decreasing the lower boundary until the function straddles
                        value = GSL_FN_EVAL( &F, low );     // get the accompanying value

                        if( fabs( low ) <= t_guess * 1e-6 || fabs( value - value_prev ) < CUTOFF )
                        {
                                std::cerr << "Couldn't adjust low. F(" << low <<
                                        ") = " << value << std::endl;
                                return low;
                        }
                        value_prev = value;
                }
                while ( value >= 0.0 );
        }


	// find the intersection on the y-axis between the random number and the function
	const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent ); // define a new solver type brent
	gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );  // make a new solver instance
									   // incl typecast?
	const Real t( findRoot( F, solver, low, high, 1e-18, 1e-12,
		"FirstPassageGreensFunction1DRad::drawTime" ) );

	return t;				// return the drawn time
}

double FirstPassageGreensFunction1DRad::drawR_f (double z, void *p)
{	struct drawR_params *params = (struct drawR_params *)p;// casts p naar type 'struct drawR_params *'
	double sum = 0;
	double An, S_Cn_An, b_An;
	int    terms = params->terms;

	for (int n=0; n<terms; n++)	// number of terms used
	{	S_Cn_An = params->S_Cn_An[n];
		b_An    = params->b_An[n];
		An      = params->An[n];
		sum += S_Cn_An * (sin(An*z) -b_An * cos(An*z) + b_An);
	}
	return sum - params->rnd;		// het snijpunt vinden met het random getal
}

const Real FirstPassageGreensFunction1DRad::drawR (const Real rnd, const Real r0, const Real t) const
{
	const Real a(this->geta());
	const Real D(this->getD());

        THROW_UNLESS( std::invalid_argument, -a <= r0 && r0 <= a );
        THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
        THROW_UNLESS( std::invalid_argument, t >= 0.0 );

	if (t == 0.0 || D == 0)
	{	return r0;	// the trivial case
	}

	struct drawR_params parameters;	// the structure to store the numbers to calculate the numbers for 1-S
	double An = 0;
	double S_Cn_An;
	double tmp0, tmp1;
	const Real L(a*2.0);
	const Real h(this->getk()/D);
	const Real S = 2.0/p_survival(t, r0);
	const Real r0_2(r0 + a);


	// produce the coefficients and the terms in the exponent and put them in the params structure
	for (int n=0; n<TERMEN; n++)
	{
		An = a_n (n+1);			// get the n-th root of tan(alfa*L)=alfa/-k
		tmp0 = An * An;			// An^2
		tmp1 = An * r0_2;			// An * z'
		S_Cn_An = S * exp(-D*tmp0*t) * (An*cos(tmp1) + h*sin(tmp1)) / (L*(tmp0 + h*h) + h);

		parameters.An[n]     = An;	// store the coefficients in the structure
		parameters.S_Cn_An[n]= S_Cn_An;	// also store the values for the exponent
		parameters.b_An[n]   = h/An;

	}
	parameters.rnd = rnd;			// store the random number for the probability
	parameters.terms = TERMEN;		// store the number of terms used

	// find the intersection on the y-axis between the random number and the function
	gsl_function F;
	F.function = &FirstPassageGreensFunction1DRad::drawR_f;
	F.params = &parameters;

	const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent ); // define a new solver type brent
	gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );  // make a new solver instance
									   // incl typecast?
	const Real z( findRoot( F, solver, 0.0, L, 1e-18, 1e-12,
		"FirstPassageGreensFunction1DRad::drawR" ) );

	return z-a;				// return the drawn place
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
