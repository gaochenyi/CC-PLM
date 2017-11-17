/**
 * Copyright (c) 2017 Chen-Yi Gao
 * 
 * LICENSE
 * ===
 * The MIT License
 * 
 * 
 * MATLAB syntax:
 * ===
 * [obj, r_grad_h_and_J] = g_r_mex_v2(...
 *   S, ...
 *   N, B, q, ...
 *   w, B_eff, ...
 *   r, r_h_and_J, ...
 *   lambda, ...
 *   SkipCheckFlag)
 * 
 *  S    uint8     [0, 255], N rows, B columns (column-major)
 *  N    uint64    length of sequence (for check: to be robust)
 *  B    uint64    number of sequences (for check: to be robust)
 *  q    uint64    $q \le 256$ since S is uint8.
 *  w    double    $\{ w_b \}$
 *  r    uint64    node index, [1,N]
 *  B_eff      double    $B_{\text{eff}} = \sum_{b=1}^B w_b$
 *  r_h_and_J  double    $q + q^2(N-1)$ rows, $1$ columns
 *  lambda     double    2 elements: first is $\lambda_h$, second is $\lambda_J$
 *   
 *  `SkipCheckFlag` is just a placeholder. 7 inputs will trigger a series of
 *  validation on inputs, while 8 inputs will skip most of the validation for
 *  speed (thus the user should be responsible for the correctness of inputs).
 *  By measuring, skipping check reduces time by $0.4\%$.
 *  
 */


#include "mex.h"
#include "g_r.v02.h"

void mexFunction(
  int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
{
  if (nlhs > 2) {
    mexErrMsgIdAndTxt(
      "g_r_mex:nlhs",
      "This function can only calculate objective and gradient.");
  }
  // default to check, 10-th arugment is a placeholder to skip check
	if (nrhs != 9 && nrhs != 10) {
    mexErrMsgIdAndTxt(
      "g_r_mex:nrhs",
      "Number of arguments supported: 9, 10\n"
      "provided: %d", nrhs);
  }

  const mxArray *pm_S           = prhs[0];
  const mxArray *pm_N           = prhs[1];
  const mxArray *pm_B           = prhs[2];
  const mxArray *pm_q           = prhs[3];
  const mxArray *pm_w           = prhs[4];
  const mxArray *pm_B_eff       = prhs[5];
  const mxArray *pm_r           = prhs[6];
  const mxArray *pm_h_r_and_J_r = prhs[7];
  const mxArray *pm_lambda      = prhs[8];

  // check when only 9 arguments are provided
  if (nrhs == 9) {
    /* type check */

    // class of mxArray
    if (   !mxIsUint8(pm_S)
        || !mxIsUint64(pm_N)
        || !mxIsUint64(pm_B)
        || !mxIsUint64(pm_q)
        || !mxIsUint64(pm_r)
        || !mxIsDouble(pm_w)
        || !mxIsDouble(pm_B_eff)
        || !mxIsDouble(pm_h_r_and_J_r)
        || !mxIsDouble(pm_lambda) )
    {
      mexErrMsgIdAndTxt(
        "g_r_mex:prhs:WrongType",
        "Requirement:\n"
        "   uint8:    S\n"
        "  uint64:    N,  B,  q,  r\n"
        "  double:    w,  B_eff,  h_r_and_J_r,  lambda");
    }

    // all inputs should be real
    if (   mxIsComplex(prhs[0])
        || mxIsComplex(prhs[1])
        || mxIsComplex(prhs[2])
        || mxIsComplex(prhs[3])
        || mxIsComplex(prhs[4])
        || mxIsComplex(prhs[5])
        || mxIsComplex(prhs[6])
        || mxIsComplex(prhs[7])
        || mxIsComplex(prhs[8]) )
    {
      mexErrMsgIdAndTxt(
        "g_r_mex:prhs:IsComplex",
        "\tAll inputs should be real.");
    }


    /* check dimensions of matrices and range of value for scalars */

    // dimension check for scalars
    if ( mxGetNumberOfElements(pm_N) != 1
      || mxGetNumberOfElements(pm_B) != 1
      || mxGetNumberOfElements(pm_q) != 1
      || mxGetNumberOfElements(pm_r) != 1
      || mxGetNumberOfElements(pm_B_eff) != 1 )
    {
      mexErrMsgIdAndTxt(
        "g_r_mex:prhs:ScalarWrong",
        "\tAll of {N, B, q, r, B_eff} should be scalar.");
    }

    // check if dimensionos match
    const size_t N0 = *((uint64_t *) mxGetData(pm_N));
    const size_t B0 = *((uint64_t *) mxGetData(pm_B));
    const size_t q0 = *((uint64_t *) mxGetData(pm_q));
    const size_t r0 = *((uint64_t *) mxGetData(pm_r));

    // S, dimension
    if ( mxGetM(pm_S) != N0
      || mxGetN(pm_S) != B0 )
    {
      mexErrMsgIdAndTxt(
        "g_r_mex:prhs:S",
        "\tColumns of `S` are considered as sequences.\n"
        "\tThus `S` should be a matrix consists of `N` rows and `B` columns");
    }

    // q, range
    if (q0 < 2 || q0 > 256) {
        mexErrMsgIdAndTxt(
          "g_r_mex:prhs:q",
          "Requirement on q:\n"
          "\tq >= 2 since 1-state Potts model is trivial.\n"
          "\tq <= 256 since `S` can only stores 0~255.\n");
    }

    // w, dimension
    if (mxGetNumberOfElements(pm_w) != B0) {
      mexErrMsgIdAndTxt(
        "g_r_mex:prhs:w",
        "\t`w`, which contains weights of sequences, mismatches `S`.\n");
    }

    // r, range
    if (r0 > N0 || r0 == 0) {
      mexErrMsgIdAndTxt(
        "g_r_mex:prhs:r",
        "\t`r` should be integers in [1,N].");
    }

    // h_r_and_J_r, dimension
    if ( mxGetM(pm_h_r_and_J_r) != (q0 + q0*q0*(N0-1))
      || mxGetN(pm_h_r_and_J_r) != 1 )
    {
      mexErrMsgIdAndTxt(
        "g_r_mex:prhs:h_r_and_J_r",
        "\t`h_r_and_J_r`, "
        "which contains h_r and J_r in pseudo-likelihood of node r, "
        "should be a (q + q*q*(N-1)) * 1 matrix\n");
    }

    // lambda, dimension and range
    if (mxGetNumberOfElements(pm_lambda) != 2) {
      mexErrMsgIdAndTxt(
        "g_r_mex:prhs:lambda",
        "\t`lambda`, which contains lambda_h and lambda_J, "
        "should contains 2 real numbers.");
    }
    const double lambda_h = mxGetPr(pm_lambda)[0];
    const double lambda_J = mxGetPr(pm_lambda)[1];
    if (lambda_h < 0.0 || lambda_J < 0.0) {
      mexErrMsgIdAndTxt("g_r_mex:prhs:lambda",
        "\t`lambda` may not be negative.");
    }
  }

  const size_t   N     = mxGetM(pm_S);
  const size_t   B     = mxGetN(pm_S);
  const size_t   q     = *((uint64_t *) mxGetData(pm_q));
  const size_t   r     = *((uint64_t *) mxGetData(pm_r));
  const uint8_t *S     = (uint8_t *) mxGetData(pm_S);
  const double  *w     = mxGetPr(pm_w);
  const double   B_eff = mxGetPr(pm_B_eff)[0];
  const double  *h_r   = mxGetPr(pm_h_r_and_J_r);
  const double  *J_r   = h_r + q;
  const double   l_h   = mxGetPr(pm_lambda)[0];
  const double   l_J   = mxGetPr(pm_lambda)[1];

  /**
   * `mxCreateDoubleMatrix` initializes each element to 0, which is required by
   * `g_r`
   */
  plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
  double *obj = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(q + q*q*(N-1), 1, mxREAL);
  double *grad_h_r = mxGetPr(plhs[1]);
  double *grad_J_r = grad_h_r + q;

  g_r(obj, grad_h_r, grad_J_r,
      B, N, q, S, w, B_eff, r-1, h_r, J_r, l_h, l_J);

}