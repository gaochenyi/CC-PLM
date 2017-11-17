% Copyright (c) 2017 Chen-Yi Gao
% 
% LICENSE
% ===
% The MIT License
% 
% DESCRIPTION
% ===
% For q-state Potts model, when $q>2$, the coupling between node $i$ and $j$ is
% parameterized by a matrices $J_{ij} (s_i, s_j)$. To get a scalar score for
% each position pair $(i,j)$, this function uses the $l_2$-norm,
% $\lVert \mathbf{x} \rVert = \sqrt{\sum_i x_i^2}$. Given a matrices describing
% coupling between node $i$ and $j$, the score is defined as
% $$
% S_{ij} = \sqrt{\sum_{s_i=2}^q \sum_{s_j=2}^q J^2_{ij}(s_i, s_j)}
% $$
% where $s_i = 1$ denotes the gap/unknown state (such as letter N in nucleic
% acid case, gap in amino/nucleic acid case).
% 
% Thus, this function assumes that first row and first column of J_ij concern
% the gap state.
% 
% OUTPUT
% ===
% `table_i_j_score` --- Every column, formatted as `[i; j; S_ij]`, contains the
% pair $(i,j)$ and score $S_{ij}$

function table_i_j_score = score_coupling_L2_no_gap(h_and_J,q,N)

if ~isa(q, 'double')
  q = double(q);
end
if ~isa(N, 'double')
  N = double(N);
end

%% Extract J_ij and calculate the score
table_i_j_score = zeros(3, N*(N-1)/2);

idx_no_gap = reshape(1:q*q, [q q]);
idx_no_gap = idx_no_gap(2:q,2:q);
idx_no_gap = idx_no_gap(:);

l = 1;
for i = 1:N-1
  J_i = reshape( h_and_J(q+1:q+q*q*(N-1), i), [q q N-1] );
  for j = i+1:N
    J_ij = J_i(:,:,j-1);  % since j > i
    % J_ij = reshape( h_and_J(q+q*q*(j-2)+1:q+q*q*(j-1),i), [q q]);
    J_ji = reshape( h_and_J(q+q*q*(i-1)+1:q+q*q*(i),j), [q q]);
    
    J_ij_sym = (J_ij + J_ji.') / 2;   % J_ji(s_j, s_i) -> J_ji(s_i, s_j)    
    score_excluding_gap = norm(J_ij_sym(idx_no_gap));
    
    table_i_j_score(:,l) = [i; j; score_excluding_gap];
    l = l+1;
  end
end

end
