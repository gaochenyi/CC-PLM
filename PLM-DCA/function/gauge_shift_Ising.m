% Copyright (c) 2017 Chen-Yi Gao
% 
% LICENSE
% ===
% The MIT License
% 
% DESCRIPTION
% ===
% Shift `h_r` and `J_r` to Ising gauge (zero-sum gauge).

function r_h_and_J_Ising = gauge_shift_Ising(r_h_and_J, q, N)

if ~isa(q, 'double')
  q = double(q);
end
if ~isa(N, 'double')
  N = double(N);
end

h_r = r_h_and_J(1:q);
J_r = reshape( r_h_and_J(q+1:end), [q q N-1] );

%% shift to Ising gauge
h_r = h_r - mean(h_r);

J_r_avg_col = mean(J_r,1);

% 2016b and later: implicit expansion
J_r = J_r - J_r_avg_col - mean(J_r,2) + mean(J_r_avg_col, 2);

% % 2016a and before
% J_r = ...
%   bsxfun(@plus, ...
%     bsxfun(@minus, ...
%       bsxfun(@minus, J_r, J_r_avg_col), ...
%     mean(J_r,2)), ...
%   mean(J_r_avg_col, 2));

r_h_and_J_Ising = [h_r; J_r(:)];

end
