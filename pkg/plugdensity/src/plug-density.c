#include <R.h>
#include <Rinternals.h>
// for M_SQRT_PI ..
#include <Rmath.h>

#include "plugdensity.h"

double plugin_h(double *x, R_xlen_t nx)
{

/************************************************************************
 *	 Version: 2022, extracted from plugin_dens() by Martin Maechler
 *
 *	 Purpose:
 *
 *	 Compute Eva Herrmann's iterative plug-in bandwidth
 *
 *	 This version only uses the gauss kernel.
 *
 *  INPUT:
 *	x[]  double    _sorted_ data array
 *	nx   int       length of  x
 *  OUTPUT:
 *	    {double}   estimated iterative plugin bandwidth
 **************************************************************************/

    const double I_7 = 1./7.,
	rt2pi = sqrt(2 * M_PI),
	rtpi2 = 2. * M_SQRT_PI,
	co2 = 1./rtpi2;

    const int iter = 5;

    /* initializations */

    double xiqr=x[(3 * nx)/4 - 1] - x[nx / 4];/* = IQR(x[]) */

    /* estimate inflation constant c */
    double n = (double) nx,
	n2 = n*n,
	h2 = (0.920 * xiqr) / pow(n, I_7),
	h3 = (0.912 * xiqr) / pow(n, 1./9.);

    double s2,s3,d2,d3, rhat2,rhat3, co1, a,t;

    s2 = s3 = 0.;
    for (R_xlen_t i = 0; i <= nx-2 ; i++) {
	for (R_xlen_t j=i+1 ; j <= nx-1 ; j++) {
	    t = x[i] - x[j];
	    a = t / h2;	    d2= a*a;
	    a = t / h3;	    d3= a*a;
	    if(d2 > 50 && d3 > 60) break;

	    s2 += exp(-d2/2.)*(  3. + d2*(-6.+d2));
	    s3 += exp(-d3/2.)*(-15. + d3*(45.+d3*(-15.+d3)));
	}
    }
    rhat2 = 2.*s2/(rt2pi*n2*pow(h2,5)) +  3./(rt2pi*n*pow(h2,5));
    rhat3 =-2.*s3/(rt2pi*n2*pow(h3,7)) + 15./(rt2pi*n*pow(h3,7));
    co1= 1.357 * pow(rhat2/rhat3, I_7);
    a = 1.132795764/(pow(rhat3, I_7)* sqrt(n));

/* MM: FIXME?  below we drop all terms  exp(-t^2 / 2) as soon as |t| > 5;
       -----   where exp(- 25/2) is "only" 3.727e-6
 */

    /* loop over iterations */

    for (int it = 1; it <= iter; it++) { // improve a
	s2=0.;
	for (R_xlen_t i = 0; i <= nx-2; i++) {
	    for (R_xlen_t j=i+1; j <= nx-1; j++) {
		t = (x[i] - x[j])/a;
		d2= t*t;
		if (d2 > 50) break; // speedup (see also above)
		s2 += exp(-d2/2.)*(3.+d2*(-6.+d2));
	    }
	}
	a = rt2pi * n * pow(a, 5);
	rhat2= 2.*s2/(n * a) + 3./a;

	/* estimate bandwidth by asymptotic formula */
	t = co2/(rhat2*n);
	a = co1*pow(t, I_7);
    }
    return /* h = */ pow(t, 0.2);
}

// FIXME: Need to use SEXP ... and .Call(...) from R  so can use large vectors !!
void plugin_dens(double *x, int *n, double *z, int *m, double *f, double *h)
{

/************************************************************************
 *	 Version: 1995
 *
 *	 Purpose:
 *
 *	 Simple	 Subroutine for kernel density estimation
 *	 with iterative plug-in bandwidth selection
 *
 *	 This version only uses the gauss kernel and estimates only
 *	 the density itself and not its derivatives.
 *
 *  INPUT:
 *	x[]  double    sorted data array
 *	n    int       length of  x
 *	z[]  double    output grid (sorted array)
 *	m    int       length of z
 *  OUTPUT:
 *	f[]  double    estimated density (array of length m)
 *	h    double    estimated iterative plugin bandwidth
 *
 *
 **************************************************************************/

    R_xlen_t nx = (R_xlen_t)*n;
    *h = plugin_h(x, nx);

    /* Estimate density f[i] = f_h(z[i]) with plugin bandwidth h : */

    const double rt2pi = sqrt(2 * M_PI);
    int jbegin = 0, jend = 0;
    jbegin = jend = 0;
    for (int i = 0; i < *m; i++) {
	double s = 0., t;
	int j;
	for(j = jbegin; j <= jend ; j++) {
	    t=(z[i] - x[j])/(*h);
	    if(t > 5. && j < nx-1) {
		jbegin++;
		continue;
	    }
	    s += exp(-t*t/2.);
	}
	for (jend = j; jend <= nx-1 ; jend++) {
	    t=(z[i] - x[jend])/(*h);
	    if(t < -5.) break;
	    s += exp(-t*t/2.);
	}
	f[i] = s/((double)nx*(*h)*rt2pi);
	jend--;
    }
}


// called from R:  .Call(h_pluginEH, x)
SEXP h_pluginEH(SEXP x_) {
    x_ = isReal(x_) ? Rf_duplicate(x_) : coerceVector(x_, REALSXP);
    PROTECT(x_);
    R_xlen_t n = XLENGTH(x_);
    R_qsort(REAL(x_), 1, (size_t)n);
    double h = plugin_h(REAL(x_), n);
    UNPROTECT(1);
    return ScalarReal(h);
}

