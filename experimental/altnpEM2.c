#include <R.h>
#include <Rmath.h>

void altnpEM2(
    int *nngrid, /* size of grid */
    int *nn, /* sample size */
    int *mm, /* number of components */
    int *rr, /* number of repeated measurements */
    double *hh, /* bandwidth */
    double *data,  /* n by r vector of observations */
    double *grid, /* grid points */
    double *old_f, /* ngrid by m by r array of current density values on grid */
    double *new_f, /* ngrid by m by r array to hold new values on grid */
    double *lambda, /* current vector of mixing parameters */
    double *post, /* n by m matrix of posterior probabilities */
    double *conv /* n by m by r working area of convolution values */
    ) {
  int n=*nn, m=*mm, r=*rr, ngrid=*nngrid, i, j, k, a, c, minj=0;
  double sum, xik, *fjk, *cvalue, two_h_squared =2*(*hh)*(*hh);
  double Delta = (grid[2]-grid[1]);
  /*  First step:  Fill in 'conv' table of convolution values. 
      This speeds up the second step a lot. */
  for (j=0; j<m; j++) {
    for (k=0; k<r; k++) {
      fjk = old_f + ngrid*(j + m*k);
      for (i=0; i<n; i++) {
        xik=data[i + n*k];
        cvalue = conv + i + n*(j + m*k);
        *cvalue = 0.0;
        for (c = 0; c<ngrid; c++) {
          *cvalue += exp(-(xik - grid[c])*(xik-grid[c])/two_h_squared) * log(fjk[c]);
        }
        *cvalue *= Delta;
      }
    }
  }
  /* Second step:  Calculate the new_f values.  Note that each of
     the ngrid * m * r entries requires a separate loop of size n.  */
  for (a=0; a<ngrid; a++) {
    cvalue = conv;
    for (k=0; k<r; k++) {
      for (j=0; j<m; j++) {
        sum = 0.0;
        for (i=0; i<n; i++) {
          xik=data[i + n*k];          
          sum += post[i + n*j] * exp(-(xik - grid[a])*(xik-grid[a])/two_h_squared) / 
          *(cvalue++);
        }
        new_f[a + ngrid*(j + m*k)] *= sum / lambda[j] / n;
      }
    }
  }
}

void altnpEM2_Estep(
    int *nngrid, /* size of grid */
    int *nn, /* sample size */
    int *mm, /* number of components */
    int *rr, /* number of repeated measurements */
    double *hh, /* bandwidth */
    double *data,  /* n by r vector of observations */
    double *grid, /* grid points */
    double *f, /* ngrid by m by r array of density values on grid */
    double *lambda, /* current vector of mixing parameters */
    double *post, /* n by m matrix of posterior probabilities */      
    double *loglik /* scalar value of penalized loglikelihood */
    ) {
  int n=*nn, m=*mm, r=*rr, ngrid=*nngrid, i, j, k, a, c, minj=0;
  double sum, conv, xik, *fjk, two_h_squared =2*(*hh)*(*hh);
  double Delta = (grid[2]-grid[1]) / *hh / sqrt(2*3.14159265358979);
  *loglik=0.0;
  for (i=0; i<n; i++) {
    sum=0.0;
    for (j=0; j<m; j++) {
      post[i + n*j] = lambda[j];
      for (k=0; k<r; k++) {
        fjk = f + ngrid*(j + m*k);
        xik=data[i + n*k];
        conv = 0.0;
        for (c = 0; c<ngrid; c++) {
          conv += exp(-(xik - grid[c])*(xik-grid[c])/two_h_squared) * fjk[c];
        }
        conv *= Delta;
        post[i + n*j] *= conv;
      }
      sum += post[i + n*j];
    }
    *loglik += log(sum);
    for(j=0; j<m; j++) {
      post[i + n*j] /= sum;
    }
  }
}


