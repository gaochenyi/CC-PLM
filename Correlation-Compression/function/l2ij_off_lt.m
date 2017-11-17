% 
% e.g. N = 3
% 
%   | 1 2 3 |
% --|-------|
% 1 | 1     |
% 2 | 2 4   |
% 3 | 3 5 6 |
% 
% l -> (i,j)
% e.g.
% 5 -> (3,2)
% 
% cf. p87 of Notebook

function ij = l2ij_off_lt(l, N)

lp = N*(N+1)/2 - l + 1;
jp = round(sqrt(2*lp));
ip = lp - jp*(jp-1)/2;
ij = [N+1-ip, N+1-jp];

end
