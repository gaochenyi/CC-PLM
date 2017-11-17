% DESCRIPTION
% ===
% This function calculates the mutual information, in unit of bit (shannon), of
% two random variables X and Y, which is denoted as I(X,Y).
%
% Assumption
% ===
% `px`, `py` are stored as column vectors while `pxy` is stored as a matrix.
% $$\sum_y p(x,y) = p(x)$$
% $$\sum_x p(x,y) = p(y)$$
% $$\sum_x \sum_y p(x,y) = 1$$
%
% **No check on input data!**
%
% FORMAT
% ===
% - px as p(x)
% - py as p(y)
% - pxy as p(x,y)
% - square_q as the total number of entries in pxy
%
% |   -   |   p(y1)   |   p(y2)   | ...  |
% | :---: | :-------: | :-------: | :--: |
% | p(x1) | p(x1, y1) | p(x1, y2) |      |
% | p(x2) | p(x2, y1) | p(x2, y2) |      |
% |  ...  |           |           |      |
%
% HISTORY
% ===
% - 2017-10-20  v3
%   - 3/2/1.5 times speed of v2a when q=2/3/5. (When q is large, v2a is faster.)
%   - interface changed: (px, py, pxy, square_q) --> (px, py, pxy, q)
%
% - 2017-10-15  v2a
%   - faster (~2x of v1)
%   - interface changed: (px, py, pxy) --> (px, py, pxy, square_q)
%
% - 2017-07-07  v1
%   - initial draft

% TEST
% ===
%     pxy = 1/4*ones(2);         % I(x,y) = 0
%     pxy = [1/2, 0; 0, 1/2];    % I(x,y) = 1
%     px = sum(pxy,2);
%     py = sum(pxy,1).';

function I = calc_MI(px, py, pxy, q)

%% computational routine
I = 0;
for i = 1:q
  for j = 1:q
    if (pxy(i,j) > 0)
      I = I + pxy(i,j) * log2( pxy(i,j) / px(i) / py(j) );
    end
  end
end

end
