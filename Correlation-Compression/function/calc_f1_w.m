% Assumption:
% ===
% `datai` is used as a 1D array. The q possible states on site i are encoded as
% integers in [1,q]. It is possible that `max(data) < q`.
%
% **No check on input data!**
%
% FORMAT
% ===
% `fi(k)` contains $f_i(k)$, the 1-point frequency on site i.
% `datai` contains $s_i^b$.
% `weights` contains weights of samples.
% $$
% f_i(k) = \frac{1}{B_eff} \sum_{b=1}^{B} w_b \delta(s_i^b,k)
% $$
% where $B_eff = \sum_b w_b$.
%
% HISTORY
% ===
% - 2017-10-15  v2
%   - faster
%   - interface changed:  (datai, weights, q) --> (datai, B, q, wiehgts, B_eff)
% 
% - 2017-10-13  v1
%   - copied from `calc_fi.m` and add support for weighted data

function fi = calc_f1_w(datai, B, q, weights, B_eff)

%% computational routine
fi = zeros(q,1);
for b = 1:B
  fi(datai(b)) = fi(datai(b)) + weights(b);
end
fi = fi/B_eff;

end
