#include <R.h>
#include <Rmath.h>

void kernelsmooth(int *nn, int *mm, int *rr, double *x, double *hh, double *z, 
  double *f, double *fk) {
  int n=*nn, m=*mm, r=*rr, rn=r*n, mn=m*n, i, jn, kn, ii;
  double sum, tmp, h=*hh, const1 = -1.0/(2.0*h*h);
  double const2 = 0.39894228040143267794/h; /* 1/(h*sqrt(2*pi)) */
 
  /* Calculate regular kernel density estimate */
  for(jn=0; jn<mn; jn+=n) {
    for(i=0; i<n; i++) {
      for(kn=0; kn<rn; kn+=n) {
        sum = 0.0;
        for(ii=0; ii<n; ii++) {
          tmp = x[i + kn]-x[ii + kn];
          sum += z[ii + jn]*exp(tmp*tmp*const1); /* Use normal kernel */
        }
        f[i + jn + kn*m] = sum * const2; 
      }
      fk[i + jn] = f[i + jn];  /* Create fkernel matrix from f array */
      for(kn=n; kn<rn; kn+=n) {
        fk[i + jn] *= f[i + jn + kn*m];
      }
    }
  }
}

void spEM_RM_general (int *n, int *m, int *r, double *x, double *h, double *z,
                      double *lam, double *A) {
  int i, j, k, a, iter, index, nn=*n, mm=*m, rr=*r, mr=mm*rr;
  double sum, diff, minus1over2h2 = -1/2/(*h)/(*h);
  double *tmpptr;

  for(iter=0; iter<1; iter++) {
    
    /* First step:  calculate marginal proportions (lambda, unnormalized) */
    for (j=0; j<mm; j++) {
      lam[j]=0.0;
      tmpptr = z + j*mm;
      for (i=0; i<nn; i++) {
        lam[j] += tmpptr[i]; /* tmpptr[i] == z[i + j*mm] */
      }
    }
    
    /* Second step:  Estimate densities (A matrix, unnormalized) */    
    for (k=0; k<rr; k++) {
      tmpptr = x + k*rr;
      for (i=0; i<nn ; i++) {
        for (j=0; j<mm; j++) {
          index = i + j*mm + k*mr;
          A[index] = 0.0;
          for (a=0; a<nn; a++) {
            diff = tmpptr[i] - tmpptr[a];  /* tmpptr[i] == x[i + k*rr] */
            A[index] += z[i + j*mm] * exp(minus1over2h2 * diff * diff);
          }
        }
      }
    }
    
    /* Third step:  calculate posteriors (z matrix, unnormalized) */    
    for (j=0; j<mm; j++) {
      tmpptr = z + j*mm;
      for (i=0; i<nn; i++) {
        tmpptr[i] = lam[j]; /* tmpptr[i] == z[i + j*mm] */
        for (k=0; k<rr; k++) {
          tmpptr[i] *= A[i + j*mm + k*mr];
        }
      }
    }

    /* Fourth step:  Normalize each row of z matrix */
    for (i=0; i<nn; i++) {
      sum = 0.0;
      for (j=0; j<mm; j++) {
        sum += z[i + j*mm];
      }
      if (sum <= 0.0) Rprintf("Warning!  Underflow error\n");
      for (j=0; j<mm; j++) {
        z[i + j*mm] /= sum;
      }
    }
  }
}

