#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdio.h>

void kernelsmoothrepeated(int *nn, int *mm, int *rr, double *x, double *hh, double *z, double *f){
  int n=*nn, i, ii;
  int rn = *rr*n, mn=*mm*n, jn, kn, kkn;
  double sum1, sum2, tmp, h=*hh, xik;
  double const1 = -0.5 / (h * h);
  double const2 = 0.39894228040143267794/(h*(double)(*rr)); /* .3989...=1/(sqrt(2*pi)) */

  for(jn=0; jn<mn; jn+=n) {
    for(i=0; i<n; i++) {
      f[i + jn] = 1.0;
      for(kn=0; kn<rn; kn+=n) {
        sum1 = 0.0;
        xik = x[i + kn];
        for(ii=0; ii<n; ii++) {
          sum2 = 0.0;
          for(kkn=0; kkn < rn; kkn+=n) {
            tmp = xik - x[ii + kkn];
            sum2 += exp(tmp * tmp * const1); /* Using normal kernel */
          }
          sum1 += z[ii + jn] * sum2;
        }
        f[i + jn] *= sum1 * const2;
      }
    }
  }
}



