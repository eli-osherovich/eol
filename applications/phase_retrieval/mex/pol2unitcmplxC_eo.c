#include <math.h>
#include "mex.h"
#include "ipp.h"


/*
  This MEX file implements MATLAB's 'pol2cart' command with a unit
  complex output, i.e., the magnitude is assumed to be unit using
  Intel's IPP library.  */


/*
  To compile run
  mex -O pol2unitcmplxC_eo.c -I/opt/intel/ipp/include -lippvm
*/


void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  double *re, *im, *magnitude, *phase;
  mwSize ndim;
  const mwSize *dims;
  size_t nelements;

  /* check proper input and output */
  if( nrhs != 1 || !mxIsDouble(prhs[0]))
    mexErrMsgTxt("One input is required.");
  
  
  nelements = mxGetNumberOfElements(prhs[0]);
  ndim = mxGetNumberOfDimensions(prhs[0]);
  dims = mxGetDimensions(prhs[0]);

  /* should test that the inputs are conforomal */
  /* meanwhile skipped */
  
  /* create output array */
  plhs[0] = mxCreateNumericArray(ndim, dims, mxDOUBLE_CLASS, mxCOMPLEX);


  re = mxGetPr(plhs[0]);
  im = mxGetPi(plhs[0]);
  phase = mxGetPr(prhs[0]);

  /* calculate cart2pol with IPP */
  ippsCos_64f_A53 (phase, re, nelements);
  ippsSin_64f_A53 (phase, im, nelements);
}
