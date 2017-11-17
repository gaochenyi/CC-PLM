% Assumption
% ===
% `datai` and `dataj` are used as 1D arrays. The q possible states on site i are
% encoded as integers in [1,q]. It is possible that `max(data) < q`.
% 
% **No check on input data!**
% 
% FORMAT
% ===
% fij(k,l) contains $f_{ij}(k,l)$, the 2-point frequency of site i and site j.
% `datai` contains $s_i^b$ while `dataj` contains $s_j^b$.
% `weights` contains weights of samples.
% $$
% f_{ij}(k,l) = \frac{1}{B_eff} \sum_{b=1}^{B} w_b \delta(s_i^b,k) \delta(s_j^b,l)
% $$
% where $B_eff = \sum_b w_b$.
% 
%     | --------------- | --------------- |
%     | f(s_i=1, s_j=1) | f(s_i=1, s_j=2) |
%     | --------------- | --------------- |
%     | f(s_i=2, s_j=1) |                 |
%     | --------------- | --------------- |
% 
% HISTORY
% ===
% - 2017-10-15  v2
%   - adapted from `calc_f2_v2`
%     - faster
%     - interface changed:
%         (datai, dataj, weights, q) --> (datai, dataj, B, q, weights, B_eff)
%     - function name changed: calc_fij_w --> calc_f2_w
% 
% - 2017-10-13  v1
%   - copied from `calc_fij.m` and add support for weighted data

function fij = calc_f2_w(datai, dataj, B, q, weights, B_eff)

%% computational routine
fij = zeros(q);
for b = 1:B
  fij(datai(b),dataj(b)) = fij(datai(b),dataj(b)) + weights(b);
end
fij = fij/B_eff;

end
