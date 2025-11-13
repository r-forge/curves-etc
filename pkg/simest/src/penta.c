#include <R_ext/Arith.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <float.h>

void penta(int dim[], double a[], double b[], double c[], double d[], 
	double e[], double g[], double u[], double v[], double w[], double p[], double q[],
	double r[], double s[], double t[], double aa[], double dd[], double z[]){
	int n = dim[0];
	float a11 = 0.0, a12 = 0.0, a21= 0.0, a22 = 0.0, b1 = 0.0, b2 = 0.0;
	for(int i = 0; i < n; i++){
		dd[i] = d[i] + e[i]*v[i];
		aa[i] = a[i] + dd[i]*v[i+1] + e[i]*w[i];
		u[i+2] = (g[i] - dd[i]*u[i+1] - e[i]*u[i])/aa[i];
		v[i+2] = -(b[i] + dd[i]*w[i+1])/aa[i];
		w[i+2] = -c[i]/aa[i];
	}
	p[1] = 1.0;
	q[0] = 1.0;
	for(int i = 0; i < n; i++){
		p[i+2] = -(dd[i]*p[i+1] + e[i]*p[i])/aa[i];
		q[i+2] = -(dd[i]*q[i+1] + e[i]*q[i])/aa[i];
	}
	r[n-2] = 1;
	s[n-1] = 1;
	for(int i = n-3; i >= 0; i--){
		t[i] = v[i+2]*t[i+1] + w[i+2]*t[i+2] + u[i+2];
		s[i] = v[i+2]*s[i+1] + w[i+2]*s[i+2] + p[i+2];
		r[i] = v[i+2]*r[i+1] + w[i+2]*r[i+2] + q[i+2];
	}
	a11 = 1 - q[n] - w[n]*r[0];
	a12 = -(p[n] + v[n] + w[n]*s[0]);
	a21 = -(v[n+1]*r[0] + w[n+1]*r[1] + q[n+1]);
	a22 = 1 - p[n+1] - v[n+1]*s[0] - w[n+1]*s[1];
	b1 = w[n]*t[0] + u[n];
	b2 = v[n+1]*t[0] + w[n+1]*t[1] + u[n+1];
	z[n-2] = (b1*a22 - b2*a12)/(a11*a22 - a12*a21);
	z[n-1] = (-b1*a21 + b2*a11)/(a11*a22 - a12*a21);
	for(int i = 0; i < n-2; i++)
		z[i] = t[i] + s[i]*z[n-1] + r[i]*z[n-2];
}
