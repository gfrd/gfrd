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
#include "Defs.hpp"

#define TERMEN 500
/*
enum EventType
{
    REACTION = 0,
    ESCAPE = 1
};
*/
class FirstPassageGreensFunction1DRad
{

public:
	FirstPassageGreensFunction1DRad(const Real D, const Real k)
		: D(D), k(k)
	{	;	// do nothing
	}

	~FirstPassageGreensFunction1DRad()
	{	;	// empty
	}

	void seta(const Real a)
	{	THROW_UNLESS( std::invalid_argument, a >= 0.0 );
		this->a =a;
	}

	const Real geta() const
	{	return this->a;
	}

	const Real getk() const
	{	return this->k;
	}

	const Real getD() const
	{	return this->D;
	}

	// Berekent de kans dat het deeltje zich nog in het domein bevindt op tijdstip t,
	// de survival probability
	const Real p_survival (const Real r0, const Real t) const;

	// Berekent de kans om het deeltje op plaats z te vinden op tijdstip t
	const Real prob_r (const Real r0, const Real r, const Real t) const;

	// Berekent de kans om het deeltje op plaats z te vinden op tijdstip t, gegeven dat het deeltje
	// zich nog in het domein bevindt.
	const Real calcpcum (const Real r0, const Real r, const Real t) const;

	// Berekent de totale kansflux die het systeem verlaat op tijdstip t
	const Real flux_tot (const Real r0, const Real t) const;

	// Berekent de kansflux die het systeem verlaat door de radiating boundary condition op x=0
	// op tijdstip t
	const Real flux_rad (const Real r0, const Real t) const;

	// Berekent de flux door de radiating boundary condition als fractie van de totale flux
	// Hiermee kunnen we de relatieve kans bepalen dat het deeltje het systeem heeft verlaten door
	// de radiating boundary en niet door de absorbing boundary.
	const Real fluxRatioRadTot (const Real r0, const Real t) const;

	const EventType drawEventType( const Real rnd, const Real r0, const Real t ) const;

	// Trekt een tijd uit de propensity function, een first passage time.
	const Real drawTime (const Real r0, const Real rnd) const;

	const Real drawR (const Real r0, const Real t, const Real rnd) const;

private:
	// Berekent de wortels van tan(aL)=-ak/h
	const Real a_n(int n) const;

	const Real An (const Real a_n, const Real r0) const;

	const Real Bn (const Real a_n, const Real r) const;

	const Real Cn (const Real a_n, const Real t) const;

	struct tan_f_params
	{	Real L;
		Real h;
	};

	static double tan_f (double x, void *p);
	// this is the appropriate definition of the function in gsl

	struct drawT_params
	{	double exponent[TERMEN];	// use 10 terms in the summation for now
		double Xn[TERMEN];
		int    terms;
		double rnd;			// the random number associated with the time
	};

	static double drawT_f (double t, void *p);

	struct drawR_params
	{	double An[TERMEN];
		double S_Cn_An[TERMEN];
		double b_An[TERMEN];
		int terms;
		double rnd;			// the random number associated with the time
	};

	static double drawR_f (double z, void *p);

	const Real D;	// The diffusion constant
	const Real k;	// The reaction constant
	Real a;		// The distance from the internal origin to one of the boundaries
			// The total domain is 2*a long
};
