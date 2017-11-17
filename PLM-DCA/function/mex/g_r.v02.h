/**
 * Copyright (c) 2017 Chen-Yi Gao
 * 
 * # License
 * 
 * The MIT License
 * 
 * 
 * # Description
 *  
 * This computational routine calculates weighted and L2 regularized objective
 * and its gradient, corresponding to asymmetric pseudo-likelihood maximisation
 * of q-state Potts model. (cf. `PLM-Potts-implement.md` or p.81 of notebook)
 * 
 * Minimization of this function corresponds to maximization of pseudo-
 * likelihood
 * 
 * 
 * # Output (by pointer)
 * 
 *  *Initialization* needed for *all* elements in obj[], grad_h_r[] and
 *  grad_J_r[]. They should be initialized to 0 before passed into this
 *  function.
 *  
 *  obj      : pointer for storage of objective value (1-D, 1 elements)
 *  grad_h_r : pointer for storage (1-D, q elements)
 *  grad_J_r : pointer for storage (1-D, q*q*(N-1) elements)
 *  
 *  
 * # Input
 * 
 *  B      : number of sequences in MSA data
 *  N      : number of nodes
 *  q      : number of states, should <= 256 since S[0] is uint8_t
 *  S      : MSA data; uint8_t, [0, q-1]; data is stored sequence-by-sequence,
 *           rather than node-by-node (1-D, N*B elements)
 *  w      : weights of sequences (1-D, B elements)
 *  B_eff  : effective number of sequences, B_eff = \sum_b w_b
 *  r      : node index, 0-indexing, [0, N-1]
 *  h_r    : local field on node r, 1-D, h_r[0]corresponds to $h_r(1)$
 *  J_r    : coupling matrix $J_{r i} \forall i \in \partial r$; matrices are
 *           stored in column-major (1-D, q*q*(N-1) elements)
 *  l_h    : lambda for L2 regularization on h_r
 *  l_J    : lambda for L2 regularization on J_r
 *  
 *  
 *  # Note for implementation
 *  
 *  In this function, `q` is chosen to be `size_t` to avoid integer overflow
 *  when calculating index from $s_i^b$, which is represented as `uint8_t`.
 *  
 *  
 *  # History
 *  
 *  ## 2017-08-09  v2
 *  - (q*S_i^b + q*q*i) -> (q*(S_i^b + q*i))
 *  - change the order of `mxFree`
 *  - minor reformatting
 *  
 *  ## 2017-08-06  v1
 */

#include <math.h>   // exp() and log()
#include <stdint.h> // uint8_t
#include "mex.h"    // mxMalloc & mxFree, mexErrMsgIdAndTxt

extern inline
void g_r(
  double *obj, double *grad_h_r, double *grad_J_r,
  const size_t B, const size_t N, const size_t q,
  const uint8_t *S,
  const double *w, const double B_eff,
  const size_t r,
  const double *h_r, const double *J_r,
  const double l_h, const double l_J)
{
  /**
   * Given a sample $\underline{s}^b$:
   * 
   * - `Num` is the numerator in marginal probability; k-th element stands for
   *   $\exp( h_r(k) + \sum_{i \neq r} J_{r i}(k, s_i^b) )$. Later this array
   *   will be used to calculate $Z_r(\underline{s}_{\partial r}^b)$.
   * 
   * - `Z_r` is denominator in marginal probability; `Z_r` equals sum of `Num`
   * 
   * - `Prob`, namely, $P(s_r^b = k | \underline{s}_{\partial r}^b)$, is the
   *   marginal probability of node $r$ being $k$. Later this array will be
   *   used to calculate gradient of h and J.
   *   
   */
  double *Num  = (double*) mxMalloc(sizeof *Num * q);
  double *Prob = (double*) mxMalloc(sizeof *Prob * q);

  // loop over samples
  for (size_t b = 0; b < B; b++) {

    /* begin: calculate $h_r(k) + \sum_{i \neq r} J_{r i}(k, s_i^b)$ */
    for (size_t k = 0; k < q; k++) {
      Num[k] = h_r[k];
    }
    // i < r
    for (size_t i = 0; i < r; i++) {
      for (size_t k = 0; k < q; k++) {
        // Num[k] += J_r[k + q*S[N*b+i] + q*q*i];
        Num[k] += J_r[k + q*(S[N*b+i] + q*i)];        // J_{ri}(k,s_i^b), i < r
      }
    }
    // i > r
    for (size_t i = r+1; i < N; i++) {
      for (size_t k = 0; k < q; k++) {
        Num[k] += J_r[k + q*(S[N*b+i] + q*(i-1))];    // J_{ri}(k,s_i^b), i > r
      }
    }
    /* end */

    /**
     * - Z_r
     * - exp(x), x = $h_r(k) + \sum_{i \neq r} J_{r i}(k, s_i^b)$
     */
    double Z_r = 0;
    for (size_t k = 0; k < q; k++) {
      if (Num[k] > 709.0) {
        mexErrMsgIdAndTxt(
          "g_rC:overflow:exp",
          "r = %ld:  "
          "`Num[k]` is too large and likely to overflow `exp(Num[k])`\n",
          r);
      }
      
      // Now Num[k] is safe
      Num[k] = exp(Num[k]);
      Z_r   += Num[k];
    }

    // now `Prob` can be got from `Num` and `Z_r`
    for (size_t k = 0; k < q; k++) {
      Prob[k] = Num[k] / Z_r;
    }

    size_t srb = S[N*b+r];  // s_r^b

    // value of objective
    *obj -= w[b] * log(Prob[srb]);

    // graddient of h_r(k)
    for (size_t k = 0; k < q; k++) {
      grad_h_r[k] += w[b] * Prob[k];
    }
    grad_h_r[srb] -= w[b];

    // graddient of J_{r j}(k, s_j^b)
    for (size_t j = 0; j < r; j++) {    // j < r
      for (size_t k = 0; k < q; k++) {
        grad_J_r[k + q*(S[N*b+j] + q*j)] += w[b]*Prob[k];
      }
      grad_J_r[srb + q*(S[N*b+j] + q*j)] -= w[b];
    }
    for (size_t j = r+1; j < N; j++) {  // j > r
      for (size_t k = 0; k < q; k++) {
        grad_J_r[k + q*(S[N*b+j] + q*(j-1))] += w[b]*Prob[k];
      }
      grad_J_r[srb + q*(S[N*b+j] + q*(j-1))] -= w[b];
    }

  }

  mxFree(Prob);
  mxFree(Num);


  /*************************************************
   *    from sum to mean (easy to be vectorized)
   *************************************************/
  *obj /= B_eff;
  for (size_t k = 0; k < q; k++) {
    grad_h_r[k] /= B_eff;
  }
  for (size_t i = 0; i < (q*q*(N-1)); i++) {
    grad_J_r[i] /= B_eff;
  }


  /*******************************************
   *    add L2 regulator
   *******************************************/
  if (l_h > 0.0) {
    for (size_t k = 0; k < q; k++) {
      grad_h_r[k] += l_h*h_r[k]*2;
      *obj        += l_h*h_r[k]*h_r[k];
    }
  }

  if (l_J > 0.0) {
    for (size_t i = 0; i < (q*q*(N-1)); i++) {
      grad_J_r[i] += l_J*J_r[i]*2;
      *obj        += l_J*J_r[i]*J_r[i];
    }
  }
}
