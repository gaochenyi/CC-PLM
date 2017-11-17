/**
 * 
 * MATLAB syntax: fij = calc_f2_w(datai, dataj, B, q, weights, B_eff)
 * 
 * # INPUT type
 * 
 *  uint8:    datai,  dataj
 * uint64:    B,  q
 * double:    w,  B_eff
 * 
 * 
 * # Assumption:
 * 
 * `datai` and `dataj` are used as 1D arrays. The q possible states are encoded
 * as integers in [1,q]. It is possible that `max(data) < q`.
 * 
 * No check on input data!!!
 * 
 * Also, this implementation assumes that size_t is equal to uint64_t.
 * 
 *
 * # FORMAT
 * 
 * fij(k,l) contains $f_{ij}(k,l)$, the 2-point frequency of site i and site j.
 * `datai` contains $s_i^b$ while `dataj` contains $s_j^b$.
 *
 * $$
 * f_{ij}(k,l) = \frac{1}{B_eff} \sum_{b=1}^{B} w_b \delta(s_i^b,k) \delta(s_j^b,l)
 * $$
 *
 * where $B_eff = \sum_b w_b$.
 *
 * | --------------- | --------------- |
 * | f(s_i=1, s_j=1) | f(s_i=1, s_j=2) |
 * | --------------- | --------------- |
 * | f(s_i=2, s_j=1) |                 |
 * | --------------- | --------------- |
 * 
 * 
 * # HISTORY
 * 
 * 2017-10-20  v1
 * 
 */

#include <cstdint>
#include "mex.h"
#include "calc_f2_w_col.hpp"

void mexFunction(
  int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
{
  if (nlhs > 1) {
    mexErrMsgIdAndTxt(
      "calc_f2_w_mex_uint8:nlhs",
      "This function produces 1 output.");
  }
  if (nrhs != 6 && nrhs != 7) { // 7-th aguments as Flag to skip check
    mexErrMsgIdAndTxt(
      "calc_f2_w_mex_uint8:nrhs",
      "Number of arguments needed: 6\n"
      "provided: %d", nrhs);
  }

  const mxArray *pm_datai = prhs[0];
  const mxArray *pm_dataj = prhs[1];
  const mxArray *pm_B     = prhs[2];
  const mxArray *pm_q     = prhs[3];
  const mxArray *pm_w     = prhs[4];
  const mxArray *pm_B_eff = prhs[5];

  if (nrhs == 6) {  // default to check
    // class of mxArray
    if (   !mxIsUint8(pm_datai)
        || !mxIsUint8(pm_dataj)
        || !mxIsUint64(pm_B)
        || !mxIsUint64(pm_q)
        || !mxIsDouble(pm_w)
        || !mxIsDouble(pm_B_eff) )
    {
      mexErrMsgIdAndTxt(
        "calc_f2_w_mex_uint8:prhs:WrongType",
        "Requirement:\n"
        "   uint8:    datai,  dataj\n"
        "  uint64:    B,  q\n"
        "  double:    w,  B_eff");
    }

    // all inputs should be real
    if (   mxIsComplex(prhs[0])
        || mxIsComplex(prhs[1])
        || mxIsComplex(prhs[2])
        || mxIsComplex(prhs[3])
        || mxIsComplex(prhs[4])
        || mxIsComplex(prhs[5]) )
    {
      mexErrMsgIdAndTxt(
        "calc_f2_w_mex_uint8:prhs:IsComplex",
        "All inputs should be real.");
    }
  }

  const size_t B = *((uint64_t *) mxGetData(pm_B));
  const size_t q = *((uint64_t *) mxGetData(pm_q));

  if (nrhs == 6) {  // default to check
    const size_t B_i = mxGetNumberOfElements(pm_datai);
    const size_t B_j = mxGetNumberOfElements(pm_dataj);
    const size_t B_w = mxGetNumberOfElements(pm_w);

    // dimension
    if (B_i != B || B_j != B) {
      mexErrMsgIdAndTxt(
        "calc_f2_w_mex_uint8:prhs:data",
        "`datai` and/or `dataj` does not match B.");
    }
    if (B_w != B) {
      mexErrMsgIdAndTxt(
        "calc_f2_w_mex_uint8:prhs:weights",
        "`weights` does not match data.");
    }

    // q, range
    if (q < 2 || q > 256) {
      mexErrMsgIdAndTxt(
        "calc_f2_w_mex_uint8:prhs:q",
        "Requirement on q:\n"
        "\tq >= 2 since 1-state case is trivial.\n"
        "\tq <= 256 since data can only stores 0~255.\n");
    }
  }

  mxArray* pm_fij = mxCreateDoubleMatrix(q, q, mxREAL);
  if (pm_fij == NULL) {
    mexErrMsgIdAndTxt(
      "calc_f2_w_mex_uint8:OutOfMemory",
      "`mxCreateDoubleMatrix` failed.");
  }
  const uint8_t* datai = (uint8_t*) mxGetData(pm_datai);
  const uint8_t* dataj = (uint8_t*) mxGetData(pm_dataj);
  const double* w = mxGetPr(pm_w);
  const double B_eff = mxGetPr(pm_B_eff)[0];
  double* fij = mxGetPr(pm_fij);
  calc_f2_w_col<uint8_t, size_t>(datai, dataj, q, B, w, B_eff, fij);
  plhs[0] = pm_fij;
}
