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


#include "FirstPassagePairGreensFunction2D3.hpp"

using namespace std;

int main (void)
{
	// initializing
	const Real D=1E-12;		// all the constants of our system
	const Real k=1E-17;	// h=k/D also sometime refered to b
	const Real sigma=1E-7;	// a
	const Real a=5E-7;	// b
	const Real r0=3E-7;

	const Real r = 4e-7;		// The radius where to draw theta's when just drawing theta;s
	Real t = 0.03;			// The time to use when drawing positions
	Real rnd = 0;
	srand(23947);                   // even de random functie initialiseren

	FirstPassagePairGreensFunction2D gf(D, k, sigma);
	gf.seta (a);


	// Producing data by drawing from the distributions
	for (int i=0; i<10000; i++)
	{
		rnd = Real(rand())/RAND_MAX;			// drawing the firstpassage times
		Real time2D = gf.drawTime (rnd, r0);
		cout << rnd << " " << time2D << " ";

		rnd = Real(rand())/RAND_MAX;			// drawing the radius for GIVEN TIME
		Real place2D = gf.drawR (rnd, r0, t);
		cout << rnd << " " << place2D << " ";

		rnd = Real(rand())/RAND_MAX;			// drawing the angle for GIVEN RADIUS
		Real angle2D = gf.drawTheta (rnd, r, r0, t);
		cout << rnd/2 << " " << angle2D << " ";

		rnd = Real(rand())/RAND_MAX;			// drawing from the 2D distribution
		     angle2D = gf.drawTheta (rnd, place2D, r0, t);
		rnd = Real(rand())/RAND_MAX;			// draw a new number to make the distribution
		if (rnd > 0.5)					// symmetric
		{	angle2D = 2*M_PI - angle2D;
		}
		cout << place2D << " " << angle2D << " ";

		rnd = Real(rand())/RAND_MAX;
		int event = gf.drawEventType (rnd, r0, t);
		cout << rnd << " " << event << " ";

		cout << endl;
	}

/*
	// Producing time dependent data
	for (t=0.0001; t < 0.3; t+= 0.001)
	{
		cout << t << " ";				// Put down the time

		Real S_t = gf.p_survival(t, r0);		// The survival probability
		cout << S_t << " ";

		Real Js_t = gf.leaves(t,r0);			// The flux through the inner boundary
		cout << Js_t << " ";

		Real Ja_t = gf.leavea(t,r0);			// The flux through the outer boundary
		cout << Ja_t << " ";

		rnd = Real(rand())/RAND_MAX;
		int event = gf.drawEventType (rnd, r0, t);	// The eventtype (escape/reaction)
		cout << event << " ";

		cout << endl;
	}
*/
}
