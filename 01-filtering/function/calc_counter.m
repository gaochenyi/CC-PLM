% Copyright (c) 2017 Chen-Yi Gao
%
% LICENSE
% ===
% See 'LICENSE.txt' in the outermost folder
%
% INPUT
% ===
% `data` contains samples of a discrete random variable X: {x_b}. X has `q`
% possible states, which is represented as integers in [1,q]. `B` is the number
% of samples. It is possible that `max(data) < q`.
%
% **No check on the range of data!**
%
% OUTPUT
% ===
% `c(k)` contains the counter for k-th state.
% $$
% c(k) = \sum_{b=1}^{B} \delta(data_b,k)
% $$
%
% HISTORY
% ===
% - 2017-10-17  v1
%   - adapted from `calc_f1.m`

function c = calc_counter(data, B, q)

%% computational routine
c = zeros(q,1);
for b = 1:B
  c(data(b)) = c(data(b)) + 1;
end

end
