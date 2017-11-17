% Copyright (c) 2017 Chen-Yi Gao
%
% LICENSE
% ===
% See 'LICENSE.txt' in the outermost folder
%
% DESCRIPTION
% ===
% minimization of g_r
%
% INPUT
% ===
% |    name    | description                                               |
% | ---------- | --------------------------------------------------------- |
% | S (uint8)  | [0,q-1], columns as sequences/samples/configurations      |
% | N (uint64) | length of sequences (number of nodes/spins)               |
% | B (uint64) | number of sequences/samples/configurations                |
% | q (uint64) | number of possible states                                 |
% | weights    | sequences can be weigted                                  |
% | B_eff      | B_eff = sum(weights)                                      |
% | r (uint64) | node index (1-based)                                      |
% | lambdas    | [lambda_h lambda_J]                                       |
% | skip       | non-zero to skip built-in check of `g_r_mex_v2`           |
% | options    | passed to `minFunc`                                       |
%
% OUTPUT
% ===
% r_h_and_J = [h_r(:); J_r(:)], where h_r and J_r are in Ising gauge.
%
% HISTORY
% ===
% - 2017-10-18  v3.0
%   - name changed from `PLM_L2_node` to `min_g_r`
%
% - 2017-10-09  v2
%   - `g_r_mex` -> `g_r_mex_v2`

function r_h_and_J = min_g_r(S,N,B,q,weights,B_eff,r,lambdas,skip,options)

if skip
  funObj = @(wr) g_r_mex_v2(S,N,B,q,weights,B_eff,r,wr,lambdas,'SkipCheckFlag');
else
  funObj = @(wr) g_r_mex_v2(S,N,B,q,weights,B_eff,r,wr,lambdas);
end

wr0 = zeros(q + q*q*(N-1), 1);
r_h_and_J = minFunc(funObj,wr0,options);

end
