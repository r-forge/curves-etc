#include <R_ext/Arith.h>
#include <R.h>
#include <math.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <float.h>

void pred(int dim[], double t[], double zhat[], double a[], double D[], double V[], double kk[]){
	int n = dim[0], r = dim[1], flag = 0;
	float trm1 = 0.0, trm2 = 0.0, trm3 = 0.0, trm4 = 0.0, trm5 = 0.0, trm6 = 0.0, trm7 = 0.0;
	float kk1 = 0.0;
	for(int j = 0; j < r; j++){
		flag = 0;
		if(kk[j] < t[0])
		{	//
			kk[j] = zhat[0] + D[0]*(kk[j] - t[0]);
			flag = 1;
		}	//
		if(flag !=1 && kk[j] >= t[0] && kk[j] <= t[1]){
			if(a[0] > 0)
				kk[j] = (a[0]*pow(kk[j] - t[0], 3)/(6*(t[1] - t[0])*(t[2] - t[0]))) +
				 D[0]*(kk[j] - t[0]) + zhat[0];
			else
				kk[j] = D[0]*(kk[j] - t[0]) + zhat[0];
			flag = 1;
		}	
		if(flag != 1 && kk[j] > t[n-2] && kk[j] <= t[n-1]){
			if(a[n-3] > 0)
				kk[j] = (a[n-3]*(kk[j] - t[n-2])*(t[n-1] - t[n-2])/(2*(t[n-1] - t[n-3]))) + 
					(a[n-3]*(pow(t[n-1] - kk[j],3) - pow(t[n-1] - t[n-2],3))/(6*(t[n-1] - t[n-2])*(t[n-1] - t[n-3]))) + 
					(D[n-2]*(kk[j] - t[n-2])) + zhat[n-2];
			else
				kk[j] = D[n-2]*(kk[j] - t[n-2]) + zhat[n-2];
			flag = 1;
		}	
		if(flag != 1 && kk[j] > t[n-1])
		{	//
			kk[j] = zhat[n-1] + D[n-1]*(kk[j] - t[n-2]);
			flag = 1;
		}	//
		for(int i = 2; i < (n-1); i++){
			if(flag !=1 && kk[j] > t[i-1] && kk[j] <= t[i]){
				if(kk[j] <= V[2*i - 1]){
           			trm1 = -a[i-2]*(pow(kk[j] - t[i], 3))/(6*(t[i] - t[i-1])*(t[i] - t[i-2]));
           			trm2 = -a[i-2]*(pow(t[i-1] - t[i], 2))/(6*(t[i] - t[i-2]));
           			trm3 = a[i-1]*(pow(kk[j] - t[i-1], 3))/(6*(t[i]-t[i-1])*(t[i+1] - t[i-1]));
           			trm4 = a[i-2]*(kk[j] - t[i-1])*(pow(t[i] - V[2*(i-1)], 2))/(2*(t[i] - t[i-1])*(t[i] - t[i-2]));
           			trm5 = -a[i-1]*(kk[j] - t[i-1])*(pow(V[2*(i-1)] - t[i-1], 2))/(2*(t[i] - t[i-1])*(t[i+1] - t[i-1]));
           			trm6 = D[i-1]*(kk[j] - t[i-1]);
           			trm7 = zhat[i-1];
           			kk[j] = trm1 + trm2 + trm3 + trm4 + trm5 + trm6 + trm7;
           			flag = 1;
           			break;
				} else{
					kk1 = V[2*i - 1];
	            	trm1 = -a[i-2]*(pow(kk1 - t[i], 3))/(6*(t[i] - t[i-1])*(t[i] - t[i-2]));
           			trm2 = -a[i-2]*(pow(t[i-1] - t[i], 2))/(6*(t[i] - t[i-2]));
           			trm3 = a[i-1]*(pow(kk1 - t[i-1], 3))/(6*(t[i]-t[i-1])*(t[i+1] - t[i-1]));
          			trm4 = a[i-2]*(kk1 - t[i-1])*(pow(t[i] - V[2*(i-1)], 2))/(2*(t[i] - t[i-1])*(t[i] - t[i-2]));
           			trm5 = -a[i-1]*(kk1 - t[i-1])*(pow(V[2*(i-1)] - t[i-1], 2))/(2*(t[i] - t[i-1])*(t[i+1] - t[i-1]));
           			trm6 = D[i-1]*(kk1 - t[i-1]);
					trm7 = zhat[i-1];
           			trm7 = trm1 + trm2 + trm3 + trm4 + trm5 + trm6 + trm7;
	           		trm1 = (a[i-2]*(pow(t[i] - V[2*(i-1)], 2) - pow(t[i] - V[2*i - 1], 2))/(2*(t[i] - t[i-1])*(t[i] - t[i-2])));
    	       		trm2 = (a[i-1]*(pow(V[2*i - 1] - t[i-1], 2) - pow(V[2*(i-1)] - t[i-1], 2))/(2*(t[i] - t[i-1])*(t[i+1] - t[i-1])));
           			kk[j] = trm7 + D[i-1]*(kk[j] - kk1) + (trm1 + trm2)*(kk[j] - V[2*i - 1]);
               		flag = 1;
               		break;	
				}
			}
		}
	}
}	
