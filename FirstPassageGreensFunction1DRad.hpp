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

#define MAX_TERMEN 100
#define MIN_TERMEN 20
#define EPSILON 1E-18

class FirstPassageGreensFunction1DRad
{

public:
	FirstPassageGreensFunction1DRad(const Real D, const Real k)
		: D(D), k(k), r0(0), L(INFINITY), factor(0)
	{	;	// do nothing
	}

	~FirstPassageGreensFunction1DRad()
	{	;	// empty
	}

	void setL(const Real L)
	{	THROW_UNLESS( std::invalid_argument, L >= 0.0 && r0 <= L);

		if ( 0.0 <= L && L < EPSILON )
		{	this->L = -1;		// just some random value to show that the domain is zero
			this->factor = 1.0;
		}
		else
		{	this->L = 1.0;
			this->factor=1.0/L;
		}
        }

	const Real getL() const
	{	return this->L;
	}

	void setr0(const Real r0)
	{	if ( this->L < 0.0 )
		{	THROW_UNLESS( std::invalid_argument, 0.0 <= r0 && r0 <= EPSILON );
			this->r0 = 0.0;
		}
		else
		{	THROW_UNLESS( std::invalid_argument, 0.0 <= r0 && r0 <= this->L );
			this->r0 = r0;
		}
	}

	const Real getr0() const
	{	return r0*factor;
	}

	const Real getk() const
	{	return this->k;
	}

	const Real getD() const
	{	return this->D*factor*factor;
	}

	// Berekent de kans dat het deeltje zich nog in het domein bevindt op tijdstip t,
	// de survival probability
	const Real p_survival (const Real t) const;

	// Berekent de kans om het deeltje op plaats z te vinden op tijdstip t, gegeven dat het deeltje
	// zich nog in het domein bevindt.
	const Real calcpcum (const Real r, const Real t) const;

	// Berekent de totale kansflux die het systeem verlaat op tijdstip t
	const Real flux_tot (const Real t) const;

	// Berekent de kansflux die het systeem verlaat door de radiating boundary condition op x=0
	// op tijdstip t
	const Real flux_rad (const Real t) const;

	// Berekent de flux door de radiating boundary condition als fractie van de totale flux
	// Hiermee kunnen we de relatieve kans bepalen dat het deeltje het systeem heeft verlaten door
	// de radiating boundary en niet door de absorbing boundary.
	const Real fluxRatioRadTot (const Real t) const;

	const EventType drawEventType( const Real rnd, const Real t ) const;

	// Trekt een tijd uit de propensity function, een first passage time.
	const Real drawTime (const Real rnd) const;

	const Real drawR (const Real rnd, const Real t) const;

private:
	// Berekent de wortels van tan(aL)=-ak/h
	const Real a_n(int n) const;

	const Real An (const Real a_n) const;

	const Real Bn (const Real a_n) const;

	const Real Cn (const Real a_n, const Real t) const;

	// Berekent de kans om het deeltje op plaats z te vinden op tijdstip t
	const Real prob_r (const Real r, const Real t) const;

	struct tan_f_params
	{	Real L;
		Real h;
	};

	static double tan_f (double x, void *p);
	// this is the appropriate definition of the function in gsl

	struct drawT_params
	{	double exponent[MAX_TERMEN];	// use 10 terms in the summation for now
		double Xn[MAX_TERMEN];
		int    terms;
		double rnd;			// the random number associated with the time
	};

	static double drawT_f (double t, void *p);

	struct drawR_params
	{	double An[MAX_TERMEN];
		double S_Cn_An[MAX_TERMEN];
		double b_An[MAX_TERMEN];
		int terms;
		double rnd;			// the random number associated with the time
	};

	static double drawR_f (double z, void *p);

	static const Real CUTOFF = 1e-10;

	const Real D;	// The diffusion constant
	const Real k;	// The reaction constant
	Real r0;
	Real L;
	Real factor;	// The scaling factor to make it work with all scales
};
