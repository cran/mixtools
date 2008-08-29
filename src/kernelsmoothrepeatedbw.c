#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdio.h>

void kernelsmoothrepeatedbw(int *nn, int *mm, int *rr, double *x, double *hh, double *z, double *f){
  int n=*nn, i, ii;
  int mn = *mm*n, rn=*rr*n, jn, kn, kkn;
  double sum1, sum2, tmp, xik;
  double const2 = 0.39894228040143267794/((double)(*rr)); /* .3989...=1/(sqrt(2*pi)) */
  double const1;

  for(jn=0; jn<mn; jn+=n, hh++) { /* at each iteration, *hh is the next bw value */
    const1 = -0.5 / (*hh * *hh);
    for(i=0; i<n; i++) {
      f[i + jn] = 1.0;
      for(kn=0; kn<rn; kn+=n) {
        sum1 = 0.0;
        xik = x[i + kn];
        for(ii=0; ii<n; ii++) {
          sum2 = 0.0;
          for(kkn=0; kkn<rn; kkn+=n) {
            tmp = xik - x[ii + kkn];
            sum2 +=  exp(tmp * tmp * const1); /* Using normal kernel */
          }
          sum1 += z[ii + jn] * sum2;
        }
        f[i + jn] *= sum1 * const2 / *hh;
      }
    }
  }
}



