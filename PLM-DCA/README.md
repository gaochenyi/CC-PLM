PLM-DCA
===
PLM-DCA is a MATLAB implementation (with companion C/C++ MEX Files) for [direct coupling analysis (DCA)][DCA] with the underlying algorithm being [pseudo-likelihood maximization (PLM)][PLM Potts].

This implementation also provides a function to perform PLM of the Potts model only: `PLM_L2_Asym`.

## License

[The MIT License](https://opensource.org/licenses/MIT)


How to use
===
  1. Call `mexAll_PLM` to compile all the MEX files.
  2. When used in 2016a and before, the implicit expansion used in `gauge_shift_Ising` should be achieved by `bsxfun`. (See the embedded comment.)
  3. Call the main function `PLM_DCA`.

## Syntax


    table_i_j_score = PLM_DCA(S,N,B,q,weights,lambda,numWorker)


### Input
  |   name    | description                              |
  | :-------: | ---------------------------------------- |
  |     S     | This $N\times B$ uint8 matrix is the MSA dataset, columns of which are sequences/samples. The $q$ possible states should be mapped to integers in $[0,q-1]$ and moreover scoring requires that the gap state should be mapped to $0$. |
  |     N     | length of sequences                      |
  |     B     | number of sequences                      |
  |     q     | number of possible states (e.g., 21 for protein problems) |
  |  weights  | weights associated with sequences/samples. If data is unweighted, just provided `ones(B,1)`. |
  |  lambda   | strength of $l_2$ regularization         |
  | numWorker | number of workers used: recommended value is the number of cores. |

### Output
`table_i_j_score`: this $3\times \frac{N(N-1)}{2}$ matrix contains the scores of all couplings. In each column, the third element is the score while the first two ones are the coupling's endpoints (1-based indexing).

### More details
For more detailed description, see comments embedded in files.


Overview of implementation
===
PLM-DCA is divided into two parts:

  1. `PLM_L2_Asym` performs [the asymmetric version of $l_2$ regularized PLM][PLM asym] for Potts model and transforms the final Potts parameters to Ising gauge;
  2. `score_coupling_L2_no_gap` scores couplings by a modified Frobenius norm which excludes contribution from the gap state. Here the requirement that gap state being mapped to $0$ is introduced. One coud provide a new function to use a different scoring scheme.

Friends of `PLM_DCA`
---

  1. `PLM_L2_Asym` and `PLM_L2_Asym_file` infer Potts parameters only; the latter makes inference for big Potts model less painful.
  2. `PLM_DCA_file` makes DCA on very large systems less painful.
  3. `example_PLM_DCA` and `example_PLM_L2_Asym` provide templates for calling provided functions.

Some references
---
  1. For a detailed description of PLM for Potts model—including regularization, Ising gauge, reweighting and scoring by Frobenius norm—one can check [Ekeberg et al. Phys. Rev. E 87, 012707 (2013).][PLM Potts]
  2. Here the [asymmetric variant][PLM asym] is used (i.e., the constraint $J_{ij} = J_{ji}$ is relaxed; then all $g_r$'s are minimized independently; afterwards ${\left( J_{ij} + J_{ji} \right)} / {2}$ is taken as $J^{\rm PLM}_{ij}$.)
  3. The gauge chosen is the Ising gauge, which aims to "put as much as possible of the Hamiltonian into local fields, and as few as possible into interaction energies" (quoted from [Weigt et al. PNAS 106, 67 (2009).][Weigt-2009]).
  4. The scoring scheme is the modified Frobenius norm, which excludes the contribution involving gap state.

[DCA]: https://en.wikipedia.org/wiki/Direct_coupling_analysis
[PLM Potts]: http://dx.doi.org/10.1103/PhysRevE.87.012707
[PLM asym]: https://doi.org/10.1016/j.jcp.2014.07.024
[Weigt-2009]: https://doi.org/10.1073/pnas.0805923106


Credits
===
  - This implementation is basically a re-implementation of [plmDCA][] and inspired a lot by it---many thanks to Magnus Ekeberg.
  - The minimization is accomplished by [minFunc][]---many thanks to Mark Schmidt.

[plmDCA]: https://github.com/magnusekeberg/plmDCA
[minFunc]: https://www.cs.ubc.ca/~schmidtm/Software/minFunc.html


Limitations
===
  1. This implementation can only handle a Potts model with $q \le 256$. This is essentially enough for DCA applications.
  2. MEX files in `minFunc` uses the old (MATLAB Version 7.2) array-handling API, which limits arrays to $2^{31}-1$ elements. As a result, without any modification, the limit of applicable systems is $\left[ q + q^2*(N-1) \right] \cdot \mathtt{Corr} \le 2^{31}-1$, where $\mathtt{Corr}$ (default to 100) is the number of corrections stored for L-BFGS method and $N$ is the number of nodes in the system. (Note that we denote it by $L$ in our paper.) To break through the limitation, one has two choices:
    - Disable MEX files by `options.useMEx = false;`. This pushes up the bound from $(2^{31}-1)$ to $2^{64}$, at the expense of more runtime. (It takes 15%  more time for a test dataset of size 81506x3145.)
    - Modify MEX files to use the new array-handling API.
