% Copyright (c) 2017 Chen-Yi Gao
%
% LICENSE
% ===
% See 'LICENSE.txt' in the outermost folder
%
% DESCRIPTION
% ===
% This function implements the 'Correlation Compression (CC)' procedure
% introduced in [a paper by C.-Y. Gao, H.-J. Zhou, and E. Aurell][link]
%
% [link]: https://doi.org/10.1103/PhysRevE.98.032407
%
% INPUT
% ===
% Following variables whose types are not specified will be checked by MATLAB.
%
% | name       | type    | note                                     |
% | ---------- | ------- | ---------------------------------------- |
% | MSA        | uint8   | [1,q], rows as sequences/samples         |
% | MSA_id     | char    | a character vector, identifier of MSA    |
% | len_seq    |         | length of sequences                      |
% | num_seq    |         | number of sequences                      |
% | q          |         | number of possible states on each locus  |
% | weights    | double  | weights of sequences/samples             |
% | num_MI     |         | loci are selected by `num_MI` largest MI |
% | numWorker  |         | number of workers in parfor              |
% | outputPath |         | path for output file storing MI table    |
% | NoLoad     | logical | true to re-calculate                     |
%
% OUTPUT
% ===
% `MSA_cc` contains loci selected from `MSA` by correlation compression.
% `idx_cc` is the indices of selected loci in `MSA` (1-based).
%
% HISTORY
% ===
% - 2017-10-24  v1

function [MSA_cc,idx_cc] = CC_MSA(MSA, MSA_id, len_seq, num_seq, q, weights, ...
  num_MI, numWorker, outputPath, NoLoad)

%% check with little overhead
% MSA
if ~isa(MSA, 'uint8')
  error('MSA should be provided as uint8. Use `MSA = uint8(MSA)` to convert.')
end
if size(MSA,2) ~= len_seq || size(MSA,1) ~= num_seq
  error('MSA should be provided as rows being sequences.')
end

% MSA_id
if ~ischar(MSA_id)
  error('`MSA_id` should be provided as a char vector.')
end

% weights
if ~isa(weights, 'double') || ~isreal(weights)
  error('weights should be provided as real numbers.')
end
[s1,s2] = size(weights);
if s1*s2 ~= num_seq || (s1 ~= 1 && s2 ~= 1)
  error('Weights for samples should be provided as a 1D array.')
end

% outputPath
if exist(outputPath,'dir') ~= 7
  error('The folder `%s` does not exist.', outputPath)
end

% NoLoad
if ~islogical(NoLoad)
  error('`NoLoad` should be provided as logical.')
end


%% check with acceptable overhead
% range of MSA
if (uint64(max(MSA(:))) > uint64(q) || uint64(min(MSA(:))) < 1)
  error('q possible states in MSA should be encoded as integers in [1,q].')
end


%% search path
if exist('calc_f2_w_mex_uint8','file') ~= 3 || exist('calc_MI.m','file') ~= 2
  addpath(genpath(pwd))
end


%% real work
[B,N] = size(MSA);
B_eff = sum(weights);

filename_MI = sprintf('%s--MI.mat', MSA_id);
filename_MI_full = fullfile(outputPath, filename_MI);

if NoLoad || exist(filename_MI_full, 'file') ~= 2
  %% calculate f_i(k)
  tic

  f1 = zeros(q,N);
  if numWorker > 1
    parfor (i = 1:N, numWorker)
      f1(:,i) = calc_f1_w(MSA(:,i), B, q, weights, B_eff);
    end
  else
    for i = 1:N
      f1(:,i) = calc_f1_w(MSA(:,i), B, q, weights, B_eff);
    end
  end

  time_f1 = toc;
  fprintf('Time for calculation of 1-point frequencies: %.2f s\n', time_f1);


  %% calculate MI
  fprintf('Calculating Mutual Information ...\n')
  tic

  q0 = uint64(q);
  B0 = uint64(B);
  N_lt = N-1; % lower triangle (off-diagonal)
  num_lt = N_lt*(N_lt+1)/2;
  list_MI  = zeros(1, num_lt);
  list_sub = zeros(2, num_lt, 'uint32');
  if numWorker > 1
    parfor (l = 1:num_lt, numWorker)
      sub = l2ij_off_lt(l,N_lt);
      j = sub(1) + 1;
      i = sub(2);
      %% given (i,j) and MSA, calculate MI(i,j)
      % datai = MSA(:,i);
      % dataj = MSA(:,j);
      fij = calc_f2_w_mex_uint8(MSA(:,i), MSA(:,j), B0, q0, weights, B_eff);
      % fi = f1(:,i);
      % fj = f1(:,j);
      % I = calc_MI(f1(:,i), f1(:,j), fij);
      list_MI(l) = calc_MI(f1(:,i), f1(:,j), fij, q);
      list_sub(:,l) = [i; j];
    end
  else
    for l = 1:num_lt
      sub = l2ij_off_lt(l,N_lt);
      j = sub(1) + 1;
      i = sub(2);
      %% given (i,j) and MSA, calculate MI(i,j)
      fij = calc_f2_w_mex_uint8(MSA(:,i), MSA(:,j), B0, q0, weights, B_eff);
      list_MI(l) = calc_MI(f1(:,i), f1(:,j), fij, q);
      list_sub(:,l) = [i; j];
    end
  end

  time_MI = toc;
  fprintf('\tFinished in %.2f s\n', time_MI);


  %% sort MI table in descending order
  fprintf('Rearranging MI table in descending order ...\n')
  tic

  [list_MI_sort, sidx_MI] = sort(list_MI, 'descend');
  list_sub_sort = list_sub(:,sidx_MI);

  time_sort = toc;
  fprintf('\tFinished in %.2f s\n', time_sort);


  %% save to file
  fprintf('Saving to file ...\n')
  tic

  save(filename_MI_full, 'list_MI_sort', 'list_sub_sort', ...
    'time_MI', 'time_sort', '-v7.3')

  time_save = toc;
  fprintf('\tFinished in %.2f s\n', time_save);


  %% subscripts of the selected top of MI
  sub_MI_top = list_sub_sort(:,1:num_MI);
else
  %% load from file
  fprintf('Loading MI table (got from previous run) ...\n')
  tic

  % load(filename_MI_table, 'list_MI_sort', 'list_sub_sort', 'time_MI')
  MatFileObj = matfile(filename_MI_full);
  sub_MI_top = MatFileObj.list_sub_sort(:,1:num_MI);

  time_load = toc;
  fprintf('\tFinished in %.2f s\n', time_load);
end


%% select loci from MSA by `num_MI` largest MI
idx_cc = union(sub_MI_top, []);
MSA_cc = MSA(:,idx_cc);


end
