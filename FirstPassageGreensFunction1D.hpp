#if !defined( __FIRSTPASSAGEGREENSFUNCTION1D_HPP )
#define __FIRSTPASSAGEGREENSFUNCTION1D_HPP

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

class FirstPassageGreensFunction1D
{
private:
	static const Real L_TYPICAL = 1E-8; // measure of 'sameness' when comparing floating points numbers
	static const Real T_TYPICAL = 1E-6; // This is a typical length scale of the system, may not be true!
	static const Real EPSILON = 1E-12;  // The typical timescale of the system, may also not be true!!
	static const Real PDENS_TYPICAL = 1; // Is 1E3 a good measure for the probability density?!

	static const unsigned int MAX_TERMEN = 500;
	static const unsigned int MIN_TERMEN = 20;

public:
        FirstPassageGreensFunction1D(const Real D)
                : D(D), L(INFINITY), r0(0), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
        {       ;       // do nothing
        }

        ~FirstPassageGreensFunction1D()
        {       ;       // empty
        }

        void setL(const Real L)			// This also sets the scale
        {       THROW_UNLESS( std::invalid_argument, L >= 0.0 && r0 <= L);

		if ( L <= EPSILON * l_scale )	// Use a typical domain size to determine if we are here
						// defining a domain of size 0.
		{	this->L = -1;		// just some random value to show that the domain is zero
//			this->l_scale = 1.0;	// don't touch the scales
		}
		else
		{	this->l_scale = L;	// set the scale to the given one
			this->t_scale = (L*L)/D;// set the typical time scale (msd = sqrt(2*d*D*t) )
			this->L = 1.0;		// L = L/l_scale
		}
        }

        const Real getL() const
        {       return this->L;
        }

        const Real getD() const
        {       return this->D/(l_scale*l_scale);
        }

	void setr0(const Real r0)
        {       if ( this->L < 0.0 )		// if the domain had zero size
                {       THROW_UNLESS( std::invalid_argument, 0.0 <= r0 && r0 <= EPSILON * l_scale );
                        this->r0 = 0.0;
                }
                else				// The normal case
                {       THROW_UNLESS( std::invalid_argument, 0.0 <= r0 && r0 <= this->L * l_scale );
                        this->r0 = r0;
                }
        }

	const Real getr0() const
	{	return this->r0/l_scale;
	}

	// Trekt een tijd uit de propensity function, een first passage time.
	const Real drawTime (const Real rnd) const;

	// Berekent een positie gegeven dat het deeltje zich nog in het domein bevindt en er twee absorbing
	// boundary conditions gelden
	const Real drawR (const Real rnd, const Real t) const;

	// Calculates the amount of flux leaving the left boundary at time t
	const Real leaves(const Real t) const;

	// Calculates the amount of flux leaving the right boundary at time t
	const Real leavea(const Real t) const;

	// Determines based on the flux ratios if the particle left the left or right boundary
	const EventType drawEventType( const Real rnd, const Real t ) const;

	// Berekent de kans dat het deeltje zich nog in het domein bevindt op tijdstip t,
	// de survival probability
	const Real p_survival (const Real t) const;

	// Berekent de kans om het deeltje op plaats z te vinden op tijdstip t, gegeven dat het deeltje
	// zich nog in het domein bevindt.
	const Real calcpcum (const Real r, const Real t) const;

private:
	struct drawT_params
	{	double exponent[MAX_TERMEN];	// use 10 terms in the summation for now
		double Xn[MAX_TERMEN];
		int    terms;
		Real tscale;
		double rnd;			// the random number associated with the time
	};

	static double drawT_f (double t, void *p);

	struct drawR_params
	{	double S_Cn_An[MAX_TERMEN];
		double n_l[MAX_TERMEN];
		int terms;
		double rnd;			// the random number associated with the time
	};

	static double drawR_f (double z, void *p);

	// Berekent de kans om het deeltje op plaats x of y te vinden op tijdstip t
	const Real prob_r (const Real r, const Real t) const;

private:
        const Real D;   // The diffusion constant
        Real L;         // The length of your domain (also the l_scale, see below)
        Real r0;
	Real l_scale;	// This is the 'length scale' of your system (1e-14 or 1e6).
                        // We scale everything to 1 with this
	Real t_scale;	// This is the time scale of the system.
};
#endif // __FIRSTPASSAGEGREENSFUNCTION1D_HPP
