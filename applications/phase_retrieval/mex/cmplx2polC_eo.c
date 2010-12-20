#include <math.h>
#include "mex.h"
#include "ipp.h"


/*
  This MEX file implements MATLAB's 'cart2pol' command for a complex
  argument (instead of two real) using Intel's IPP library.  */


/*
  To compile run
  mex -O cmplx2polC_eo.c -L	\
  /opt/intel/Compiler/11.1/072/ipp/em64t/sharedlib/ -lippsem64t
*/


void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double *re, *im, *magnitude, *phase;
  mwSize ndim;
  const mwSize *dims;
  size_t nelements;

  /* check proper input and output */
  if( nrhs != 1 || !mxIsComplex(prhs[0]))
    mexErrMsgTxt("One complex input required.");
  /* else if(nlhs != 2) */
  /*   mexErrMsgTxt("Two output arguments expected."); */
  
  nelements = mxGetNumberOfElements(prhs[0]);
  ndim = mxGetNumberOfDimensions(prhs[0]);
  dims = mxGetDimensions(prhs[0]);
  
  /* create output arrays */
  plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxREAL);

  re = mxGetPr(prhs[0]);
  im = mxGetPi(prhs[0]);
  magnitude = mxGetPr(plhs[1]);
  phase = mxGetPr(plhs[0]);

  /* calculate cart2pol with IPP */
  ippsCartToPolar_64f(re, im, magnitude, phase, nelements);
}
