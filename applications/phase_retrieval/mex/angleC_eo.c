#include <math.h>
#include "mex.h"
#include "ipp.h"


/* This MEX file implements MATLAB's 'angle' command using Intel's IPP
   library.
   There are several available functions:

   ippsAtan2_64f_A50 - provides 50 correct bits
   ippsAtan2_64f_A53 - provides 53 correct bits
   ippsPhase_64f - precision is not indicated however, it seems to
   provide about 12 decimal digits */


/*
  To compile run
  mex -O angleC_eo.c -I/opt/intel/ipp/include -lippvm 
*/


void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double *re, *im, *result;
  mwSize ndim;
  const mwSize *dims;
  int nelements;

  /* check proper input and output */
  if( nrhs != 1 || !mxIsComplex(prhs[0]))
    mexErrMsgTxt("One complex input required.");
  else if(nlhs > 1)
    mexErrMsgTxt("Too many output arguments.");
  
  nelements = mxGetNumberOfElements(prhs[0]);
  ndim = mxGetNumberOfDimensions(prhs[0]);
  dims = mxGetDimensions(prhs[0]);
  
  /* create output array */
  plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);

  re = mxGetPr(prhs[0]);
  im = mxGetPi(prhs[0]);
  result = mxGetPr(plhs[0]);

  /* calculate atan2 with IPP */
  ippsAtan2_64f_A50(im, re, result, nelements);
  /* ippsPhase_64f(re, im, result, nelements); */
}
