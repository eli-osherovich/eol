#include <math.h>
#include <string.h>
#include "mex.h"
#include "blas.h"


/* See  lbfgsUpdate_eo.m for details! */

/*
  NOTE!
  
  The function does not check inputs. MATLAB will crash if wrong
  inputs are provided!

  The update is done inplace: using undocumented features of MATLAB
*/

/* To compile run
   mex -O -largeArrayDims lbfgsUpdateInplaceC_eo.c -lmwblas
*/

/* declare some MATLAB functions that are not exposed by its interface */
bool mxIsSharedArray(const mxArray *array_ptr);
mxArray *mxCreateSharedDataCopy(const mxArray *pr);
bool mxUnshareArray(const mxArray *pr, const bool noDeepCopy);

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  
  /* Variable Declarations */

  /*
    [PrevSteps, PrevGrads,H0, validIdx, wrapAround] = ...
    lbfgsUpdate_eo(validIdx, wrapAround, x_diff, grad_diff, PrevSteps, PrevGrads, H0)
  */
  

  
  double *validIdx, *wrapAround, *x_diff_r, *x_diff_i, *grad_diff_r, *grad_diff_i, *H0;

  /* PervSteps and PrevGrads are Cell arrays */
  mxArray *PrevSteps, *PrevGrads;



  /* local variables */
  double gradstep;

  mwSignedIndex   nPrev,  nVars, one = 1;
  mxArray *tmp_cell_s, *tmp_cell_g;
  double *tmp_data;
  bool x_diff_complex, grad_diff_complex;
  
  /* prevent compiler warning on unused variables */
  (void) nrhs;
  (void) plhs;

  /* minimal arguments check */
  if (0 != nlhs) {
    mexErrMsgTxt("Inplace update - no outputs required.");
  }

  /* check if altered arrays are shared */
  if (mxIsSharedArray(prhs[0])) {
    mexWarnMsgTxt("validIdx is shared");
    mxUnshareArray(prhs[0], 0);
  }
  if (mxIsSharedArray(prhs[1])) {
    mexWarnMsgTxt("wrapAround is shared");
    mxUnshareArray(prhs[1], 0);
  }
  if (mxIsSharedArray(prhs[4])) {
    mexWarnMsgTxt("PrevSteps is shared");
    mxUnshareArray(prhs[4], 0);
  }
  if (mxIsSharedArray(prhs[5])) {
    mexWarnMsgTxt("PrevGrads is shared");
    mxUnshareArray(prhs[5], 0);
  }
  if (mxIsSharedArray(prhs[6])) {
    mexWarnMsgTxt("H0 is shared");
    mxUnshareArray(prhs[6], 0);
  }
  
  


  
  /* Get Input Pointers */
  validIdx = mxGetPr(prhs[0]);
  wrapAround = mxGetPr(prhs[1]);
  x_diff_r = mxGetPr(prhs[2]);
  x_diff_i = mxGetPi(prhs[2]);
  x_diff_complex = mxIsComplex(prhs[2]);
  grad_diff_r = mxGetPr(prhs[3]);
  grad_diff_i = mxGetPi(prhs[3]);
  grad_diff_complex = mxIsComplex(prhs[3]);
  PrevSteps = prhs[4];
  PrevGrads = prhs[5];
  H0 = mxGetPr(prhs[6]);

  /* Compute number of variables (from x_diff) */
  nVars = mxGetNumberOfElements(prhs[2]);
 
  /* gradstep = real(grad_diff'*x_diff); */
  gradstep = ddot(&nVars, grad_diff_r, &one, x_diff_r, &one);
  if (x_diff_complex && grad_diff_complex) {
    gradstep += ddot(&nVars, grad_diff_i, &one, x_diff_i, &one);
  }

  /* if gradstep > 1e-10 */
  if (gradstep > 1e-10) {
    /* nPrev = length(PrevSteps); */
    nPrev = mxGetNumberOfElements(prhs[4]);
    /* if validIdx >= nPrev % it should never be greater then... */
    if (*validIdx >= nPrev) {
      /* wrapAround = 1; */
      /* validIdx = 1; */
      *wrapAround = 1;
      *validIdx = 1;
    }
    else {
      /* validIdx = validIdx + 1; */
      *validIdx += 1;
    }
    /* PrevSteps{validIdx} = x_diff; */
    /* PrevGrads{validIdx} = grad_diff; */
    tmp_cell_s = mxGetCell(PrevSteps, (mwIndex)(*validIdx) - 1);
    tmp_cell_g = mxGetCell(PrevGrads, (mwIndex)(*validIdx) - 1);

    if (NULL == tmp_cell_s) {
      /* non-populated cell. create new array */
      mxArray *tmp_array;
      if (x_diff_complex) {
	tmp_array = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxCOMPLEX);
      }
      else {
	tmp_array = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
      }
      /* allocate data with 16 byte alignment */
      tmp_data = mxCalloc(nVars*sizeof(double)/16 + 1, 16);
      memcpy(tmp_data, x_diff_r, nVars*sizeof(double));
      mxSetPr(tmp_array, tmp_data);
      if (x_diff_complex) {
	tmp_data = mxCalloc(nVars*sizeof(double)/16 + 1, 16);
	memcpy(tmp_data, x_diff_i, nVars*sizeof(double));
	mxSetPi(tmp_array, tmp_data);
      }
      mxSetM(tmp_array, nVars);
      mxSetN(tmp_array, 1);
      mxSetCell(PrevSteps, (mwIndex)(*validIdx) - 1, tmp_array);

      if (grad_diff_complex) {
	tmp_array = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxCOMPLEX);
      }
      else {
	tmp_array = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
      }
      tmp_data = mxCalloc(nVars*sizeof(double)/16 + 1, 16);
      memcpy(tmp_data, grad_diff_r, nVars*sizeof(double));
      mxSetPr(tmp_array, tmp_data);
      if (grad_diff_complex) {
	tmp_data = mxCalloc(nVars*sizeof(double)/16 + 1, 16);
	memcpy(tmp_data, grad_diff_i, nVars*sizeof(double));
	mxSetPi(tmp_array, tmp_data);
      }
      mxSetM(tmp_array, nVars);
      mxSetN(tmp_array, 1);
      mxSetCell(PrevGrads, (mwIndex)(*validIdx) - 1, tmp_array);
    }
    else {
      /* an existing cell. copy new data. */
      memcpy(mxGetPr(tmp_cell_s), x_diff_r, nVars*sizeof(double));
      if (x_diff_complex) {
	memcpy(mxGetPi(tmp_cell_s), x_diff_i, nVars*sizeof(double));
      }
      memcpy(mxGetPr(tmp_cell_g), grad_diff_r, nVars*sizeof(double));
      if (grad_diff_complex) {
	memcpy(mxGetPi(tmp_cell_g), grad_diff_i, nVars*sizeof(double));
      }
    }
      

    /* H0 = gradstep/real(grad_diff'*grad_diff); */
    if (grad_diff_complex) {
      *H0 = gradstep / (ddot(&nVars, grad_diff_r, &one, grad_diff_r, &one) +
			ddot(&nVars, grad_diff_i, &one, grad_diff_i, &one));
    }
    else {
      *H0 = gradstep / ddot(&nVars, grad_diff_r, &one, grad_diff_r, &one);
    }
  }
}
