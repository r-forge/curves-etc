#include <R_ext/Arith.h>
#include <R.h>
#include <math.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <float.h>

void derivcvxpen(int dim[], double t[], double zhat[], double a[], double D[], double V[], double kk[]){
	int n = dim[0], r = dim[1], flag = 0;
	for(int j = 0; j < r; j++){
		flag = 0;
		if(kk[j] <= t[0])
		{	//
			kk[j] = D[0];
			flag = 1;
		}	//
		if(flag !=1 && kk[j] >= t[0] && kk[j] < t[1]){
			if(a[0] > 0)
				kk[j] = (a[0]*pow(kk[j] - t[0], 2)/(2*(t[1] - t[0])*(t[2] - t[0]))) +
				 D[0];
			else
				kk[j] = D[0];
			flag = 1;
		}	
		if(flag != 1 && kk[j] >= t[n-2] && kk[j] < t[n-1]){
			if(a[n-3] > 0)
				kk[j] = (a[n-3]*(pow(t[n-1] - t[n-2],2) - pow(t[n-1] - kk[j],2))/(2*(t[n-1] - t[n-2])*(t[n-1] - t[n-3]))) + 
					D[n-2];
			else
				kk[j] = D[n-2];
			flag = 1;
		}	
		if(flag != 1 && kk[j] >= t[n-1])
		{	//
			kk[j] = D[n-1];
			flag = 1;
		}	//
		for(int i = 2; i < (n-1); i++){
			if(flag !=1 && kk[j] >= t[i-1] && kk[j] < t[i]){
				if(kk[j] <= V[2*i - 1]){
           			kk[j] = D[i-1] +
    					a[i-1]*(pow(kk[j] - t[i-1],2) - pow(V[2*(i-1)] - t[i-1],2))/(2*(t[i] - t[i-1])*(t[i+1] - t[i-1])) -
    					a[i-1]*(pow(t[i] - kk[j],2) - pow(t[i] - V[2*(i-1)],2))/(2*(t[i] - t[i-1])*(t[i] - t[i-2]));
           			flag = 1;
           			break;
				} else{
					kk[j] = D[i-1] +
		    			a[i-1]*(pow(V[2*i-1] - t[i-1],2) - pow(V[2*(i-1)] - t[i-1],2))/(2*(t[i] - t[i-1])*(t[i+1] - t[i-1])) -
    					a[i-1]*(pow(t[i] - V[2*i-1],2) - pow(t[i] - V[2*(i-1)],2))/(2*(t[i] - t[i-1])*(t[i] - t[i-2]));
               		flag = 1;
               		break;	
				}
			}
		}
	}
}	
