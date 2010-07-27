#include <math.h>
#include <string.h>
#include "mex.h"
#include "blas.h"


/* See  lbfgsDirPersist_eo.m for details! */

/*
  Please Note!

  The function does not check inputs. MATLAB will crash if wrong
  inputs are prvided!

  Due to persistent (static) memory, this function is not
  re-entrant. Moreover, you cannot run two LBFGS algorithms in
  parallel.
*/

/* To compile run
   mex -O -largeArrayDims lbfgsDirPersistC_eo.c -lmwblas
*/

static double *rho = NULL;
static mwSignedIndex nPrev_old = 0;

static void
freeMem (void)
{
  if (NULL != rho)
    {
      mxFree(rho);
    }
}

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

  /* Variable Declarations */

  /*
    d = lbfgsDir_eo(validIdx, wrapAround,...
    grad, PrevSteps, PrevGrads, H0)
  */

  mwSignedIndex validIdx;
  int wrapAround;
  double *grad_r, *grad_i, *H0;
  bool complex_flag;

  /* PervSteps and PrevGrads are Cell arrays */
  const mxArray *PrevSteps, *PrevGrads;

  /* local variables */
  double *alpha, *beta, *d_r, *d_i;

  double temp;

  mwSize lhs_dims[2];
  mwSignedIndex  i, nPrev,  nVars, one = 1;
  mxArray *tmp_cell_s, *tmp_cell_g;

  /* prevent compiler warining on unused varialbes */
  (void) nlhs;
  (void) nrhs;

  /* Get Input Pointers */

  validIdx = (mwSignedIndex) mxGetScalar(prhs[0]);
  wrapAround = 0 != mxGetScalar(prhs[1]);
  grad_r = mxGetPr(prhs[2]);
  grad_i = mxGetPi(prhs[2]);
  complex_flag = mxIsComplex(prhs[2]);
  PrevSteps = prhs[3];
  PrevGrads = prhs[4];
  H0 = mxGetPr(prhs[5]);

 
  /* Compute number of variables (from grad) and number of steps
     remembered (from PrevSteps) */
  nVars = mxGetNumberOfElements(prhs[2]);
  nPrev = mxGetNumberOfElements(prhs[3]);


  /* Allocate rho (if necessary) and make it persitent */
  if (nPrev != nPrev_old)
    {
      freeMem();
      nPrev_old = nPrev;
      rho = (double*) mxMalloc(nPrev*sizeof(double));
      mexMakeMemoryPersistent(rho);

      /* register clean-up routine */
      mexAtExit(freeMem);
    }


  /* allocate tmp Variables */
  if (wrapAround)
    {
      alpha = (double*) mxMalloc(nPrev*sizeof(double));
      beta = (double*) mxMalloc(nPrev*sizeof(double));
    }
  else
    {
      alpha =(double*) mxMalloc(validIdx*sizeof(double));
      beta = (double*) mxMalloc(validIdx*sizeof(double));
    }



  /* Intialize output vector */
  lhs_dims[0] = nVars;
  lhs_dims[1] = 1;
  if (complex_flag)
    {
      plhs[0] = mxCreateNumericArray(2,lhs_dims,mxDOUBLE_CLASS,mxCOMPLEX);
    }
  else
    {
      plhs[0] = mxCreateNumericArray(2,lhs_dims,mxDOUBLE_CLASS,mxREAL);
    }
  d_r = mxGetPr(plhs[0]);
  d_i = mxGetPi(plhs[0]);


  /* d = grad; */
  memcpy(d_r, grad_r, nVars*sizeof(double));
  if (complex_flag)
    {
      memcpy(d_i, grad_i, nVars*sizeof(double));
    }

  /* if validIdx == 0 we shall use grad as the search direction */
  if (validIdx)
    {
      /*  rho(i) = 1/real(PrevSteps{i}'*PrevGrads{i}); */
      tmp_cell_s = mxGetCell(PrevSteps, validIdx - 1);
      tmp_cell_g = mxGetCell(PrevGrads, validIdx - 1);
      temp = ddot(&nVars, mxGetPr(tmp_cell_s), &one, mxGetPr(tmp_cell_g), &one);
      if (complex_flag)
        {
          temp += ddot(&nVars, mxGetPi(tmp_cell_s), &one, mxGetPi(tmp_cell_g), &one);
        }
      rho[validIdx-1] = 1.0/temp;
    }
  for ( i = validIdx-1; i >= 0; --i )
    {
      /* alpha(i) = real(PrevSteps{i}'*d) * rho(i); */
      tmp_cell_s = mxGetCell(PrevSteps, i);
      alpha[i] = ddot(&nVars, d_r, &one, mxGetPr(tmp_cell_s), &one);
      if (complex_flag)
        {
          alpha[i] += ddot(&nVars, d_i, &one, mxGetPi(tmp_cell_s), &one);
        }
      alpha[i] *= rho[i];

      /* d = d - alpha(i) * PrevGrads{i}; */
      tmp_cell_g = mxGetCell(PrevGrads, i);
      temp = -alpha[i];
      daxpy(&nVars, &temp, mxGetPr(tmp_cell_g), &one, d_r, &one);
      if (complex_flag)
        {
          daxpy(&nVars, &temp, mxGetPi(tmp_cell_g), &one, d_i, &one);
        }
    }
  if (wrapAround)
    {
      for ( i = nPrev-1; i >= validIdx; --i )
        {
          /* alpha(i) = real(PrevSteps{i}'*d) * rho(i); */
          tmp_cell_s = mxGetCell(PrevSteps, i);
          alpha[i] = ddot(&nVars, d_r, &one, mxGetPr(tmp_cell_s), &one);
          if (complex_flag)
            {
              alpha[i] += ddot(&nVars, d_i, &one, mxGetPi(tmp_cell_s), &one);
            }
          alpha[i] *= rho[i];

          /* d = d - alpha(i) * PrevGrads{i}; */
          tmp_cell_g = mxGetCell(PrevGrads, i);
          temp = -alpha[i];
          daxpy(&nVars, &temp, mxGetPr(tmp_cell_g), &one, d_r, &one);
          if (complex_flag)
            {
              daxpy(&nVars, &temp, mxGetPi(tmp_cell_g), &one, d_i, &one);
            }
        }
    }

  /* d = H0 * d; */
  dscal(&nVars, H0, d_r, &one);
  if (complex_flag)
    {
      dscal(&nVars, H0, d_i, &one);
    }

  if ( wrapAround )
    {
      for ( i = validIdx; i < nPrev; ++i )
        {
          /* beta(i) = real(PrevGrads{i}' * d) * rho(i); */
          tmp_cell_g = mxGetCell(PrevGrads, i);
          beta[i] = ddot(&nVars, d_r, &one, mxGetPr(tmp_cell_g) , &one);
          if (complex_flag)
            {
              beta[i] += ddot(&nVars, d_i, &one, mxGetPi(tmp_cell_g) , &one);
            }
          beta[i] *= rho[i];

          /* d = d + PrevSteps{i} * (alpha(i)-beta(i)); */
          temp = alpha[i]-beta[i];
          tmp_cell_s = mxGetCell(PrevSteps, i);
          daxpy(&nVars, &temp, mxGetPr(tmp_cell_s), &one, d_r, &one);
          if (complex_flag)
            {
              daxpy(&nVars, &temp, mxGetPi(tmp_cell_s), &one, d_i, &one);
            }
        }
    }


  for ( i = 0; i < validIdx; ++i )
    {
      /* beta(i) = real(PrevGrads{i}' * d) * rho(i); */
      tmp_cell_g = mxGetCell(PrevGrads, i);
      beta[i] = ddot(&nVars, d_r, &one, mxGetPr(tmp_cell_g) , &one);
      if (complex_flag)
        {
          beta[i] += ddot(&nVars, d_i, &one, mxGetPi(tmp_cell_g) , &one);
        }
      beta[i] *= rho[i];

      /* d = d + PrevSteps{i} * (alpha(i)-beta(i)); */
      temp = alpha[i]-beta[i];
      tmp_cell_s = mxGetCell(PrevSteps, i);
      daxpy(&nVars, &temp, mxGetPr(tmp_cell_s), &one, d_r, &one);
      if (complex_flag)
        {
          daxpy(&nVars, &temp, mxGetPi(tmp_cell_s), &one, d_i, &one);
        }
    }

  /* Free Memory */
  mxFree(beta);
  mxFree(alpha);
}

