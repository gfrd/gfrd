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

class FirstPassageGreensFunction1D
{

public:
        FirstPassageGreensFunction1D(const Real D)
                : D(D)
        {       ;       // do nothing
        }

        ~FirstPassageGreensFunction1D()
        {       ;       // empty
        }

        void seta(const Real a)
        {       this->a =a;
        }

        const Real geta() const
        {       return this->a;
        }

        const Real getD() const
        {       return this->D;
        }

	// Trekt een tijd uit de propensity function, een first passage time.
	const Real drawTime (const Real rnd) const;

	// Berekent een positie gegeven dat het deeltje zich nog in het domein bevindt en er twee absorbing
	// boundary conditions gelden
	const Real drawR (const Real t, const Real rnd) const;



	// Berekent de kans dat het deeltje zich nog in het domein bevindt op tijdstip t,
	// de survival probability
	const Real p_survival (const Real t) const;

	// Berekent de kans om het deeltje op plaats x of y te vinden op tijdstip t
	const Real prob_r (const Real r, const Real t) const;

	// Berekent de kans om het deeltje op plaats z te vinden op tijdstip t, gegeven dat het deeltje
	// zich nog in het domein bevindt.
	const Real calcpcum (const Real r, const Real t) const;

private:
	struct drawT_params
	{	double exponent[TERMEN];	// use 10 terms in the summation for now
		double Xn[TERMEN];
		int    terms;
		double rnd;			// the random number associated with the time
	};

	static double drawT_f (double t, void *p);

	struct drawR_params
	{	double S_Cn_An[TERMEN];
		double n_l[TERMEN];
		int terms;
		double rnd;			// the random number associated with the time
	};

	static double drawR_f (double z, void *p);

        const Real D;   // The diffusion constant
        Real a;         // The distance from the internal origin to one of the boundaries
                        // The total domain is 2*a long
};

