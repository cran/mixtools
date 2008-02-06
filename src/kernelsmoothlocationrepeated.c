#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdio.h>

void kernelsmoothlocationrepeated(int *n, int *m, int *r, double *mu, double *x, double *h, double *z, double *f, double *f2) {
  int nn=*n, mm=*m, rr=*r, i, j, k, ii, jj, kk;
  double sum1, sum2, u1, u2, tmp1, tmp2, hh=*h, fikj;
  
  /* Calculate symmetrized f a=ii and b=kk */
  for(kk=0; kk<mm; kk++) {
    for(ii=0; ii<nn; ii++) {
      fikj = 1.0;
      for(jj=0; jj<rr; jj++){
	u1 = x[ii + jj*nn] - mu[kk];
	sum1 = 0.0;
	for(i=0; i<nn; i++) {
	  for(k=0; k<mm; k++) {
	    sum2 = 0.0;
	    for(j=0; j<rr; j++) {
	      u2 = x[i + j*nn] - mu[k];
	      tmp1 = u1 - u2;
	      tmp2 = -u1 - u2;
	      /* Use normal kernel */
	      sum2 += dnorm(tmp1/hh, 0, 1, 0) + dnorm(tmp2/hh, 0, 1, 0); /* Sum of kernel for each j=1,...,r */
	    }
	    sum1 += z[i + k*nn] * sum2;
	  }
	}
	f2[ii + jj*nn + kk*nn*rr] = sum1 / (2.0 * (double) nn * hh * (double) rr);
	fikj *= sum1 / (2.0 * (double) nn * hh * (double) rr);
      }
      f[ii + kk*nn] = fikj;
    }
  }
}
