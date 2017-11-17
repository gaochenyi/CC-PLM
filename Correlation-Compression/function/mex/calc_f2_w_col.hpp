#ifndef CALC_F2_W_COL_HPP
#define CALC_F2_W_COL_HPP

/**
 * # DESCRIPTION
 * 
 * Two discrete random variables X and Y, have q possible states. Given M
 * observations of (X,Y), which is denoted as {x_i} and {y_i}, this function
 * calculates the 2-point frequency f_{xy}. Samples are weighted---there is a
 * weight w_m associated with sample m.
 * 
 * 
 * # INPUT
 * 
 * - `datax` and `datay` contain the M samples of X and Y. q possible states are
 *   encoded as integers in [1,q].
 * - `w` contains {w_m} for M samples.
 * - M_eff = \sum_m w_m
 * - `fxy` is a double array, all elements of which have been initialized to
 *   zero. And f_{xy} is stored in *column-major* order.
 * 
 */

template<class T, class idxType>
extern inline
void calc_f2_w_col(
  const T* datax, const T* datay, const idxType q, const idxType M,
  const double* w, const double M_eff,
  double *fxy)
{
  /* fxy is initialized to all zero by user */
  // for (idxType i = 0; i < q*q; i++) { // clean fxy
  //   fxy[i] = 0.0;
  // }

  for (idxType m = 0; m < M; m++) {
    // fxy[ (i-1) + (j-1)*q ] += w[m];
    fxy[ datay[m]*q - q + datax[m] - 1 ] += w[m];
  }

  for (idxType i = 0; i < q*q; i++) {
    fxy[i] /= double(M_eff);
  }
}

#endif // CALC_F2_W_COL_HPP
