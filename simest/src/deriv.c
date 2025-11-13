#include <R_ext/Arith.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <float.h>

void deriv(int dim[], double zhat[], double ba[], double V[], double t[], double D[]){
  int n = dim[0];
  float tmp = 0.0;
  if(ba[0] > 0) 
    tmp = ba[0];
  D[0] = ((zhat[1] - zhat[0])/(t[1] - t[0])) - (tmp*(t[1] - t[0])/(6*(t[2] - t[0])));
  D[1] = D[0] + (tmp*(t[1] - t[0])/(2*(t[2] - t[0])));
  for(int i = 2; i < (n-1); i++){
    D[i] = D[i-1] 
      - ba[i-2]*(t[i] - V[(i-1)*2 + 1])*(t[i] - V[(i-1)*2 + 1])/(2*(t[i]-t[i-1])*(t[i]-t[i-2]))
      + ba[i-2]*(t[i] - V[(i-1)*2])*(t[i] - V[(i-1)*2])/(2*(t[i]-t[i-1])*(t[i]-t[i-2]))
      + ba[i-1]*(V[(i-1)*2 + 1]-t[i-1])*(V[(i-1)*2 + 1] - t[i-1])/(2*(t[i] - t[i-1])*(t[i+1]-t[i-1]))
      - ba[i-1]*(V[(i-1)*2]-t[i-1])*(V[(i-1)*2] - t[i-1])/(2*(t[i] - t[i-1])*(t[i+1]-t[i-1]));
  }
  if(ba[n-3] > 0)
    tmp = ba[n-3];
  else
    tmp = 0.0;
  D[n-1] = D[n-2] + (tmp*(t[n-1] - t[n-2])/(2*(t[n-1] - t[n-3])));
}
