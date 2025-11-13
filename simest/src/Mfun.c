#include <R_ext/Arith.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <float.h>

void Mfun(int dim[], double a[], double t[], double V[], double M[]) {
    int N = dim[0];
    int Vrow = N + 1, Vcol = 2, Mcol = N;
    float t1 = 0.0, t2 = 0.0, ts = 0.0, tri = 3.0;
    float denom1 = 1.0, denom2 = 1.0, f1num = 0.0, f2num = 0.0;
    float num = 0.0, den = 1.0;
    if(a[0] > 0){
    	V[0] = t[0];
    	V[1] = t[1];
    } else {
    	V[0] = t[0];
    	V[1] = t[0];
    }
    if(a[N-1] > 0){
    	V[(Vrow-1)*Vcol] = t[N];
    	V[(Vrow-1)*Vcol + 1] = t[N+1]; 
    } else{
    	V[(Vrow-1)*Vcol] = t[N];
    	V[(Vrow-1)*Vcol + 1] = t[N];
    }
    for(int i = 1; i < N; i++){
    	t1 = a[i]*(t[i+1] - t[i-1]);
    	t2 = a[i-1]*(t[i+2] - t[i]);
    	ts = ((t1*t[i]) - (t2*t[i+1]))/(t1 - t2);
    	if(a[i] > 0 && a[i-1] > 0){
    		V[i*Vcol] = t[i];
    		V[i*Vcol + 1] = t[i+1];
    	}
    	if(a[i] > 0 && a[i-1] < 0){
    		V[i*Vcol] = ts;
    		V[i*Vcol + 1] = t[i+1];
    	}
    	if(a[i] < 0 && a[i-1] > 0){
    		V[i*Vcol] = t[i];
    		V[i*Vcol + 1] = ts;
    	}
    	if(a[i] < 0 && a[i-1] < 0){
    		V[i*Vcol] = t[i];
    		V[i*Vcol + 1] = t[i];
    	}
    }
    for(int i = N-1; i--;){
    	denom1 = 3*(t[i+2] - t[i])*(t[i+2] - t[i])*(t[i+1] - t[i])*(t[i+1] - t[i]);
    	denom2 = 3*(t[i+2] - t[i])*(t[i+2]-t[i])*(t[i+2] - t[i+1])*(t[i+2] - t[i+1]);
    	f1num = pow(V[i*Vcol + 1] - t[i], tri) - pow(V[i*Vcol] - t[i], tri);
    	f2num = pow(V[(i+1)*Vcol + 1] - t[i+2], tri) - pow(V[(i+1)*Vcol] - t[i+2], tri);
    	M[i*Mcol + i] = (f1num/denom1) + (f2num/denom2);
    }
    for(int i = N-2; i--; ){
    	den = -(t[i+2] - t[i])*(t[i+3] - t[i+1])*(t[i+2] - t[i+1])*(t[i+2] - t[i+1]);
    	num = (pow(V[(i+1)*Vcol + 1], tri) - pow(V[(i+1)*Vcol], tri))/3
    		- (t[i+1] + t[i+2])*(pow(V[(i+1)*Vcol + 1], 2.0) - pow(V[(i+1)*Vcol], 2.0))/2
    		+ t[i+1]*t[i+2]*(V[(i+1)*Vcol + 1] - V[(i+1)*Vcol]);
    	M[i*Mcol + i + 1] = num/den;
    	M[(i+1)*Mcol + i] = M[i*Mcol + i + 1];
    }
}
