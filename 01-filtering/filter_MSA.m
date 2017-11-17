% Copyright (c) 2017 Chen-Yi Gao
%
% LICENSE
% ===
% See 'LICENSE.txt' in the outermost folder
%
% DESCRIPTION
% ===
% Implementation of filtering introduced in PLOS Genetics paper of
% [Skwark et al.](https://doi.org/10.1371/journal.pgen.1006508)
%
% INPUT
% ===
% `MSA`, a uint8 matrix, consists of `B` rows as sequences, the length of which
% is `N`. NACGT are mapped to 12345.
%
% OUTPUT
% ===
% `MSA_f` is the filtered MSA containing loci satisfying the three criteria.
% `idx` is the indices of loci selected (1-based). `MSA_f` should contains
% only N/major/minor (as 123).
%
% HISTORY
% ===
% - 2017-10-17  v1
%   - adapted from `01filter.cpp`

function [MSA_f, idx_f] = filter_MSA(...
  MSA, num_seq, len_seq, q, ...
  letter_N_max, MAF_min)

% search path
if exist('filter_locus','file') ~= 2
  addpath(genpath(pwd))
end

%% check with very little overhead
% check assumption
if q ~= 5 || ~isa(MSA, 'uint8')
  error('Assumption fails.')
end

% check the orientation of MSA
[B,N] = size(MSA);
if num_seq ~= B || len_seq ~=N
  error('MSA should be a matrix with rows being sequences.');
end


%% check with acceptable overhead
% check range of MSA
if min(MSA(:)) < 1 || max(MSA(:)) > q
  error('MSA is assumed to be encoded by integers in [1,q].')
end


%% real work
fprintf('Filtering loci ...\n');
tic

MSA_new_rep = zeros(B,N,'uint8');
labels = zeros(1,N);
parfor i = 1:N
  % see `filter_locus` for details
  [MSA_new_rep(:,i), labels(i)] = filter_locus(...
    MSA(:,i), B, q, letter_N_max, MAF_min);
end

% position of loci selected (in MATLAB-style indexing)
% 1. bi-allelic;
% 2. not so may uncertain letter: N;
% 3. not so few minor letter
idx_f = find(labels == 0).';

% filtered MSA
MSA_f = MSA_new_rep(:,idx_f);

% check with acceptable overhead
if min(MSA_f(:)) < 1 || max(MSA_f(:)) > 3
  error('`MSA_f` exceeds [1,3].')
end

% counter for each type of loci
numbers = zeros(8,1);
for i = 1:N
  numbers(labels(i)+1) = numbers(labels(i)+1) + 1;
end

time_filter = toc;

fprintf('\tFinished in %.2f s.\n', time_filter);

%
fprintf('Numbers for 8 types:\n')
for i = 0:1
  for j = 0:1
    for k = 0:1
      fprintf('%d %d %d\t%d\n',i,j,k,numbers(4*i+2*j+k+1));
    end
  end
end

end
