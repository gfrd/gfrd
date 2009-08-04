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

#define TERMEN 500
#define MINTERMEN 20
#define THRESHOLD 1E-9

// Berekent de kans dat het deeltje zich nog in het domein bevindt op tijdstip t, de survival probability
const Real FirstPassageGreensFunction1D::p_survival (const Real t) const
{
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );

	const Real a(this->geta());
	const Real D(this->getD());
	Real sum = 0, term = 0;
	Real nPI;
	Real expo(-D*t/(4.0*a*a));	// exponent -D n^2 PI^2 t / l^2
	Real n=1;

	do
	{	nPI = n*M_PI;
		term = exp(nPI*nPI*expo) * sin(nPI/2.0) * (1.0 - cos(nPI)) / nPI;
		sum += term;
		n++;
	}
	while (fabs(term/sum) > THRESHOLD || n < MINTERMEN );
		// Is this a good measure or will this fail at some point?

	return sum*2.0;
}

// Berekent de kans om het deeltje op plaats x of y te vinden op tijdstip t
const Real FirstPassageGreensFunction1D::prob_r (const Real r, const Real t) const
{
	const Real a(this->geta());
	const Real D(this->getD());

	THROW_UNLESS( std::invalid_argument, -a <= r && r <= a );
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );

	const Real expo(-D*t/(4*a*a));
	int n=1;
	Real nPI;
	Real sum = 0, term = 0;

	do
	{	nPI = n*M_PI;
		term = exp(nPI*nPI*expo) * sin(nPI/2) * sin(nPI*(r+a)/(2*a));
		sum += term;
		n +=2;		// only for odd n is sin(nPI/2) unequal to zero
	}
	while (fabs(term/sum) > THRESHOLD);

	return sum/a;
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

	const Real a(this->geta());
	const Real D(this->getD());

	const Real expo(-D/(4.0*a*a));

	struct drawT_params parameters;	// the structure to store the numbers to calculate the numbers for 1-S
	Real nPI;
	Real Xn, exponent;
	int n = 0;

	// produce the coefficients and the terms in the exponent and put them in the params structure
	do
	{
		nPI = ((Real)(n+1))*M_PI;		// 
		Xn = sin(nPI/2.0) * (1.0 - cos(nPI)) / nPI; 
		exponent = nPI*nPI*expo;

		parameters.Xn[n] = Xn;			// store the coefficients in the structure
		parameters.exponent[n]=exponent;	// also store the values for the exponent
		n++;
	}
	while (n<TERMEN);	// modify this later to include a cutoff when changes are small

	parameters.rnd = rnd;			// store the random number for the probability
	parameters.terms = TERMEN;		// store the number of terms used

/*
for (double t=0; t<0.1; t += 0.0001)
{
	std::cout << t << " " << drawT_f (t, &parameters) << std::endl;
}
*/

	// find the intersection on the y-axis between the random number and the function
	gsl_function F;
	F.function = &drawT_f;
	F.params = &parameters;

	const gsl_root_fsolver_type* solverType( gsl_root_fsolver_brent ); // define a new solver type brent
	gsl_root_fsolver* solver( gsl_root_fsolver_alloc( solverType ) );  // make a new solver instance
									   // incl typecast?
	const Real t( findRoot( F, solver, 1e-5, 1.0, 1e-18, 1e-12,
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
const Real FirstPassageGreensFunction1D::drawR (const Real t, const Real rnd) const
{
	THROW_UNLESS( std::invalid_argument, 0.0 <= rnd && rnd < 1.0 );
	THROW_UNLESS( std::invalid_argument, t >= 0.0 );

	const Real a(this->geta());
	const Real D(this->getD());

	struct drawR_params parameters;	// structure to store the numbers to calculate numbers for 1-S
	Real S_Cn_An;
	Real nPI;
	int n=0;
	const Real S = 2.0/p_survival(t);
	const Real expo = -D*t/(4.0*a*a);


//std::cout << S << std::endl;
	// produce the coefficients and the terms in the exponent and put them in the params structure
	do
	{
		nPI = ((Real)(n+1))*M_PI;			// 
		S_Cn_An = S * exp(nPI*nPI*expo) * sin(nPI/2.0) / nPI;

		parameters.S_Cn_An[n]= S_Cn_An;		// also store the values for the exponent
		parameters.n_l[n]    = nPI/(2.0*a);
		n++;
	}
	while (n<TERMEN);

	parameters.rnd = rnd ;			// store the random number for the probability
	parameters.terms = TERMEN;		// store the number of terms used

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
	const Real r( findRoot( F, solver, 0.0, 2.0*a, 1e-18, 1e-12,
		"FirstPassageGreensFunction1D::drawR" ) );

	return r-a;				// return the drawn time
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
