% Assumption
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
% $$
% f_i(k) = \frac{1}{B} \sum_{b=1}^{B} \delta(s_i^b,k)
% $$
%
% HISTORY
% ===
% - 2017-10-15  v2
%   - faster
%   - interface changed: (datai, q) --> (datai, B, q)
%   - function name changed: calc_fi --> calc_f1
% 
% - 2017-07-06  v1
%   - change interface: Encoding changes from [0,1,2] to [1,q]

function fi = calc_f1(datai, B, q)

%% computational routine
fi = zeros(q,1);
for b = 1:B
  fi(datai(b)) = fi(datai(b)) + 1;
end
fi = fi/B;

end
