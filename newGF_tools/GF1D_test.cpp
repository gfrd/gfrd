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


#include "FirstPassageGreensFunction1D.hpp"

using namespace std;

int main (void)
{
	// initializing
	const Real D=1E-12;		// all the constants of our system
	const Real a=5E-7;

	Real t = 0.03;			// The time to use when drawing positions
	Real rnd = 0;
	srand(23947);                   // even de random functie initialiseren

	FirstPassageGreensFunction1D gf(D);
	gf.seta (a);


/*
	// Producing data by drawing from the distributions
	for (int i=0; i<1000; i++)
	{
		rnd = Real(rand())/RAND_MAX;			// drawing the firstpassage times
		cout << rnd << " ";
		Real time2D = gf.drawTime (rnd);
		cout << time2D << " ";

		rnd = Real(rand())/RAND_MAX;			// drawing the radius for GIVEN TIME
		Real place2D = gf.drawR (t, rnd);
		cout << rnd << " " << place2D << " ";

		cout << endl;
	}
*/
/*
	// Producing time dependent data
	for (t=0; t < 0.3; t+= 0.001)
	{
		cout << t << " ";				// Put down the time

		Real S_t = gf.p_survival(t);		// The survival probability
		cout << S_t << " ";

		cout << endl;
	}
*/

	// producing the place dependent data
	for (Real x=-a; x<a; x+=1e-8)
	{
		cout << x << " ";

		Real p = gf.prob_r (x, t);
		cout << p << " ";

		cout << endl;
	}

}
