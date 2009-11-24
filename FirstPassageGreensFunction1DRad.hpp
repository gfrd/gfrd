#if !defined( __FIRSTPASSAGEGREENSFUNCTION1DRAD_HPP )
#define __FIRSTPASSAGEGREENSFUNCTION1DRAD_HPP

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

class FirstPassageGreensFunction1DRad
{
private:
        static const Real L_TYPICAL = 1E-8;  // This is a typical length scale of the system, may not be true!
        static const Real T_TYPICAL = 1E-6;  // The typical timescale of the system, may also not be true!!
        static const Real EPSILON = 1E-10;   // measure of 'sameness' when comparing floating points numbers
        static const Real PDENS_TYPICAL = 1; // Is 1E3 a good measure for the probability density?!

        static const int MAX_TERMEN = 500;   // The maximum number of terms used in calculating the sum
        static const int MIN_TERMEN = 20;    // The minimum number of terms

public:
	FirstPassageGreensFunction1DRad(const Real D, const Real k)
		: D(D), k(k), r0(0), L(INFINITY), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
	{	;	// do nothing
	}

	~FirstPassageGreensFunction1DRad()
	{	;	// empty
	}

	void setL(const Real L)			// This also sets the scale
	{	THROW_UNLESS( std::invalid_argument, L >= 0.0 && r0 <= L);

		if ( L < EPSILON * l_scale)	// Use a typical domain size to determine if we are here
						// defining a domain of size 0.
		{	this->L = -1;		// just some random value to show that the domain is zero
//			this->l_scale = 1.0;
		}
		else
		{	this->l_scale = L;	// set the l_scale to the given one
			this->t_scale = (L*L)/D;// set the typical time scale (msd = sqrt(2*d*D*t) )
			this->L = 1.0;		// L = L/l_scale
		}
        }

	const Real getL() const
	{	return this->L;
	}

	void setr0(const Real r0)
	{	if ( this->L < 0.0 )	// if the domain had zero size
		{	THROW_UNLESS( std::invalid_argument, 0.0 <= r0 && r0 <= EPSILON * l_scale );
			this->r0 = 0.0;
		}
		else
		{	THROW_UNLESS( std::invalid_argument, 0.0 <= r0 && r0 <= this->L * l_scale );
			this->r0 = r0;
		}
	}

	const Real getr0() const
	{	return r0/l_scale;
	}

	const Real getk() const
	{	return this->k/l_scale;		// don't forget to scale the k as well!
	}

	const Real getD() const
	{	return this->D/(l_scale*l_scale);
	}
	// Calculates the probability density of finding the particle at location z at timepoint t, given
	// that the particle is still in the domain.
	const Real calcpcum (const Real r, const Real t) const;

	// Determine which event has occured, an escape or a reaction. Based on the fluxes through the
	// boundaries at the given time. Beware: if t is not a first passage time you still get an answer!
	const EventType drawEventType( const Real rnd, const Real t ) const;

	// Draws the first passage time from the propensity function
	const Real drawTime (const Real rnd) const;

	// Draws the position of the particle at a given time, assuming that the particle is still in the
	// domain
	const Real drawR (const Real rnd, const Real t) const;


// These methods are both public and private, they are used by public methods but can also be called
// from the 'outside'. This is mainly because of debugging purposes.


	// Calculates the probability of finding the particle inside the domain at time t
	// -> the survival probability
	const Real p_survival (const Real t) const;

	// Calculates the total probability flux leaving the domain at time t
	const Real flux_tot (const Real t) const;

	// Calculates the probability flux leaving the domain through the radiative boundary at time t
	const Real flux_rad (const Real t) const;

	// Calculates the flux leaving the domain through the radiative boundary as a fraction of the
	// total flux. This is the probability that the particle left the domain through the radiative
	// boundary instead of the absorbing boundary.
	const Real fluxRatioRadTot (const Real t) const;

// End of public/private mix methods


//private:
	// Calculates the roots of tan(aL)=-ak/h
	const Real a_n(int n) const;
private:

	const Real An (const Real a_n) const;

	const Real Bn (const Real a_n) const;

	const Real Cn (const Real a_n, const Real t) const;

	// Calculates the probability density of finding the particle at location r at time t.
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
		Real   tscale;			// the timescale used for convergence
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

//	static const Real CUTOFF = 1e-10;

	const Real D;	// The diffusion constant
	const Real k;	// The reaction constant
	Real r0;
	Real L;		// The length of your domain (also the l_scale, see below)
	Real l_scale;	// This is the 'scale' of your system (1e-14 or 1e6).
			// We scale everything to 1 with this
        Real t_scale;   // This is the time scale of the system.
};
#endif // __FIRSTPASSAGEGREENSFUNCTION1DRAD_HPP
