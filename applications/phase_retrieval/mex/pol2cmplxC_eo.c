#include <math.h>
#include "mex.h"
#include "ipp.h"


/*
  This MEX file implements MATLAB's 'pol2cart' command with a complex
  output (instead of two real) using Intel's IPP library.  */


/*
  To compile run
  mex -O pol2cmplxC_eo.c -L	\
  /opt/intel/Compiler/11.1/072/ipp/em64t/sharedlib/ -lippsem64t
*/


void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double *re, *im, *magnitude, *phase;
  mwSize ndim;
  const mwSize *dims;
  size_t nelements;

  /* check proper input and output */
  if( nrhs != 2 || !mxIsDouble(prhs[0]))
    mexErrMsgTxt("Two inputs are required.");
  
  
  nelements = mxGetNumberOfElements(prhs[0]);
  ndim = mxGetNumberOfDimensions(prhs[0]);
  dims = mxGetDimensions(prhs[0]);

  /* should test that the inputs are conforomal */
  /* meanwhile skipped */
  
  /* create output array */
  plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxCOMPLEX);


  re = mxGetPr(plhs[0]);
  im = mxGetPi(plhs[0]);
  magnitude = mxGetPr(prhs[1]);
  phase = mxGetPr(prhs[0]);

  /* calculate cart2pol with IPP */
  ippsPolarToCart_64f(magnitude, phase, re, im, nelements);
}
