#include <R.h>
#include <Rmath.h>

void kernelsmoothlocation_old /* See modified version below */
(int *n, int *m, double *mu, double *x, double *h, double *z, double *f) {
  int nn=*n, mm=*m, i, j, a, b;
  double sum, u1, u2, tmp1, tmp2, hh=*h;
  
  /* Calculate symmetrized f */
  for(a=0; a<nn; a++) {
    for(b=0; b<mm; b++) {
      sum = 0.0;
      for(i=0; i<nn; i++) {
	for(j=0; j<mm; j++) {
	  u1 = (x[a]-mu[b]);
	  u2 = (x[i]-mu[j]);
	  tmp1 = u1-u2;
	  tmp2 = -u1-u2;
	  /* Use normal kernel */
	  sum += z[i + j*nn]*(dnorm(tmp1/hh, 0, 1, 0) + dnorm(tmp2/hh, 0, 1, 0));
	}
      }
      f[a + b*nn] = sum / (2.0 * (double) nn * hh);
    }
  }
}


void kernelsmoothlocation(int *n, int *m, double *mu, double *x, double *h, double *z, double *f) {
  int nn=*n, mm=*m, i, j, a, b;
  double sum, u1, u2, tmp1, tmp2, hh=*h;
  double const1 = -1.0 / (2.0 * hh * hh);
  double const2 = 0.39894228040143267794/(2.0*hh*(double)nn); /* .3989...=1/(sqrt(2*pi)) */
  
  /* Calculate symmetrized f */
  for(a=0; a<nn; a++) {
    for(b=0; b<mm; b++) {
      sum = 0.0;
      for(i=0; i<nn; i++) {
        for(j=0; j<mm; j++) {
          u1 = (x[a]-mu[b]);
          u2 = (x[i]-mu[j]);
          tmp1 = u1-u2;
          tmp2 = -u1-u2;
          /* Use normal kernel */
          sum += z[i + j*nn] * (exp(tmp1 * tmp1 * const1) + 
                                exp(tmp2 * tmp2 * const1));
        }
      }
      f[a + b*nn] = sum * const2;
    }
  }
}
