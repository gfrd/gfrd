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
    // This is a typical length scale of the system, may not be true!
    static const Real L_TYPICAL = 1E-8;
    // The typical timescale of the system, may also not be true!!
    static const Real T_TYPICAL = 1E-6;
    // measure of 'sameness' when comparing floating points numbers
    static const Real EPSILON = 1E-12;
    //E3; Is 1E3 a good measure for the probability density?!
    static const Real PDENS_TYPICAL = 1;
    // The maximum number of terms in the sum
    static const int MAX_TERMEN = 500;
    // The minimum
    static const int MIN_TERMEN = 20;

public:
    FirstPassageGreensFunction1D(const Real D)
	: D(D), L(INFINITY), r0(0), l_scale(L_TYPICAL), t_scale(T_TYPICAL)
    {
	;   // do nothing
    }

    ~FirstPassageGreensFunction1D()
    { 
	;   // empty
    }

    // This also sets the scale
    void setL(const Real L)
    {
	THROW_UNLESS( std::invalid_argument, L >= 0.0 && r0 <= L);

	// Use a typical domain size to determine if we are here
	// defining a domain of size 0.
	if ( L <= EPSILON * l_scale )
	{
	    // just some random value to show that the domain is 
	    // zero
	    this->L = -1;
	    // don't touch the scales
	    //this->l_scale = 1.0;
	}
	else
	{
	    // set the scale to the given one
	    this->l_scale = L;
	    // set the typical time scale (msd = sqrt(2*d*D*t) )
	    this->t_scale = (L*L)/D;
	    // L = L/l_scale
	    this->L = 1.0;
	}
    }

    const Real getL() const
    {
	return this->L;
    }

    const Real getD() const
    {
	return this->D/(l_scale*l_scale);
    }

    void setr0(const Real r0)
    {
	if ( this->L < 0.0 )
	{
	    // if the domain had zero size    
	    THROW_UNLESS( std::invalid_argument, 0.0 <= r0 && r0 <= EPSILON * l_scale );
	    this->r0 = 0.0;
	}
	else
	{
	    // The normal case
	    THROW_UNLESS( std::invalid_argument, 0.0 <= r0 && r0 <= this->L * l_scale );
	    this->r0 = r0;
	}
    }

    const Real getr0() const
    {
	return this->r0/l_scale;
    }

    // Draws the first passage time from the propensity function
    const Real drawTime (const Real rnd) const;

    // Draws the position of the particle at a given time, assuming that 
    // the particle is still in the
    // domain
    const Real drawR (const Real rnd, const Real t) const;

    // Calculates the amount of flux leaving the left boundary at time t
    const Real leaves(const Real t) const;

    // Calculates the amount of flux leaving the right boundary at time t
    const Real leavea(const Real t) const;

    // Determines based on the flux ratios if the particle left the left 
    // or right boundary
    const EventType drawEventType( const Real rnd, const Real t ) const;

    // Calculates the probability of finding the particle inside the 
    // domain at time t so, the survival probability
    const Real p_survival (const Real t) const;

    // Calculates the probability density of finding the particle at 
    // location z at timepoint t, given that the particle is still in the 
    // domain.
    const Real calcpcum (const Real r, const Real t) const;

private:
    struct drawT_params
    {
	// use 10 terms in the summation for now
	double exponent[MAX_TERMEN];
	double Xn[MAX_TERMEN];
	int    terms;
	Real tscale;
	// the random number associated with the time
	double rnd;
    };

    static double drawT_f (double t, void *p);

    struct drawR_params
    {
	double S_Cn_An[MAX_TERMEN];
	double n_l[MAX_TERMEN];
	int terms;
	// the random number associated with the time
	double rnd;
    };

    static double drawR_f (double z, void *p);

    // Calculates the probability density of finding the particle at 
    // location r at time t.
    const Real prob_r (const Real r, const Real t) const;

private:
    // The diffusion constant
    const Real D;
    // The length of your domain (also the l_scale, see below)
    Real L;
    Real r0;
    // This is the 'length scale' of your system (1e-14 or 1e6).
    // We scale everything to 1 with this
    Real l_scale;
    // This is the time scale of the system.   
    Real t_scale;
};
#endif // __FIRSTPASSAGEGREENSFUNCTION1D_HPP
