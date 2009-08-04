#include <sstream>
#include <iostream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_bessel.h>


#include "../FirstPassageGreensFunction1DRad.hpp"

using namespace std;

int main (void)
{
	// initializing
	const Real D=1E-12;		// all the constants of our system
	const Real k=1E-17;
	const Real a=5E-7;
	const Real r0= 0;

	Real t = 0.03;			// The time to use when drawing positions
	Real rnd = 0;
	srand(23947);                   // even de random functie initialiseren

	FirstPassageGreensFunction1DRad gf(D, k);
	gf.seta (a);


	// Producing data by drawing from the distributions
	for (int i=0; i<1000; i++)
	{
		rnd = Real(rand())/RAND_MAX;			// drawing the firstpassage times
		cout << rnd << " ";
		Real time2D = gf.drawTime (rnd, r0);
		cout << time2D << " ";

		rnd = Real(rand())/RAND_MAX;			// drawing the radius for GIVEN TIME
		cout << rnd << " ";
		Real place2D = gf.drawR (rnd, r0, t);
		cout << place2D << " ";

		rnd = Real(rand())/RAND_MAX;
		int event = gf.drawEventType (rnd, r0, t);
		cout << rnd << " " << event << " ";

		cout << endl;
	}

/*
	// Producing time dependent data
	for (t=0.001; t < 0.3; t+= 0.001)
	{
		cout << t << " ";				// Put down the time

		Real S_t = gf.p_survival(t, r0);		// The survival probability
		cout << S_t << " ";

		Real Js_t = gf.flux_rad(r0, t);			// The flux through the radiating boundary
		cout << Js_t << " ";

		Real Jtot_t = gf.flux_tot(r0, t) - Js_t;		// The total flux 
		cout << Jtot_t << " ";

		rnd = Real(rand())/RAND_MAX;
		int event = gf.drawEventType (rnd, r0, t);	// The eventtype (escape/reaction)
		cout << event << " ";

		cout << endl;
	}
*/
/*
	// producing the place dependent data
	for (Real x=-a+1E-10; x<a; x+=1e-8)
	{
		cout << x << " ";

		Real p = gf.prob_r (r0, x, t);
		cout << p << " ";

		cout << endl;
	}
*/
}
