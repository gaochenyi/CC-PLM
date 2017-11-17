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
% $$
% f_{ij}(k,l) = \frac{1}{B} \sum_{b=1}^{B} \delta(s_i^b,k) \delta(s_j^b,l)
% $$
% 
%     | --------------- | --------------- |
%     | f(s_i=1, s_j=1) | f(s_i=1, s_j=2) |
%     | --------------- | --------------- |
%     | f(s_i=2, s_j=1) |                 |
%     | --------------- | --------------- |
% 
% HISTORY
% ===
% - 2017-10-15  v2b
%   - faster
%   - function name changed: calc_fij --> calc_f2
%
% - 2017-10-15  v2a
%   - faster
%   - interface changed
%
% - 2017-07-07  v1
%   - initial draft

function fij = calc_f2(datai, dataj, B, q)

%% computational routine
fij = zeros(q);
for b = 1:B
  fij(datai(b),dataj(b)) = fij(datai(b),dataj(b)) + 1;
end
fij = fij/B;

end
