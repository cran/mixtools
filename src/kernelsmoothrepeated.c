#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdio.h>

void kernelsmoothrepeated(int *nn, int *mm, int *rr, double *x, double *hh, double *z, double *f){
  int n=*nn, m=*mm, r=*rr, i, k, j, ii, kk;
  double sum1, sum2, tmp, h=*hh, xik;
  double const1 = -1.0 / (2.0 * h * h);
  double const2 = 0.39894228040143267794/(h*(double)r); /* .3989...=1/(sqrt(2*pi)) */

  for(j=0; j<m; j++) {
    for(i=0; i<n; i++) {
      f[i + j*n] = 1.0;
      for(k=0; k<r; k++) {
        sum1 = 0.0;
        xik = x[i + k*n];
        for(ii=0; ii<n; ii++) {
          sum2 = 0.0;
          for(kk=0; kk<r; kk++) {
            tmp = xik - x[ii + kk*n];
            sum2 += exp(tmp * tmp * const1); /* Using normal kernel */
          }
          sum1 += z[ii + j*n] * sum2 * const2;
        }
        f[i + j*n] *= sum1;
      }
    }
  }
}



