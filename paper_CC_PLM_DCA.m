% HISTORY
% ===
% - 2017-10-24  v1

% MEMO
% ===
% `fastafile`   checked by `filter_FASTA`
% `dataID`      checked here
% `lambda`      checked here
% `outputPath`  checked here
% `NoLoad`      checked here

function paper_CC_PLM_DCA(fastafile, dataID, num_MI, lambda, ...
  numWorker, outputPath, NoLoad)


%% check with little overhead

% dataID
if ~ischar(dataID)
  error('`dataID` should be provided as a char vector.')
end

% lambda
if lambda < 0.0
  error('lambda should be non-negative.')
end

% outputPath
if ~ischar(outputPath)
  error('`dataPath` should be provided as a char vector.')
end
if exist(outputPath,'dir') ~= 7
  error('The folder `%s` does not exist.', outputPath)
end

% NoLoad
if ~islogical(NoLoad)
  error('`NoLoad` should be provided as logical.')
end


%% If no pool exists, a new one is created.
poolobj = gcp('nocreate');
if isempty(poolobj)
  parpool(numWorker);
end


%% search path
addpath(genpath(pwd))


%% Filtering
letter_N_max = 500;
MAF_min = 0.01;

filename_filter = sprintf(...
  '%s--gapMax_%g-MAFmin_%g-N_1-major_2-minor_3.mat', dataID, ...
  letter_N_max, MAF_min);
filename_filter_full = fullfile(outputPath, filename_filter);

% 1. Loading results of previous run is not allowed.
% 2. Results of previous do not exist.
if NoLoad || exist(filename_filter_full, 'file') ~= 2
  [MSA_f, idx_f] = filter_FASTA(fastafile, letter_N_max, MAF_min);
  save(filename_filter_full, 'MSA_f', 'idx_f', '-v7.3');
else
  load(filename_filter_full, 'MSA_f', 'idx_f');
end

% remove duplicate samples 
MSA_f_unique = unique(MSA_f,'rows');

% MSA info
[B_f,N_f] = size(MSA_f_unique);
q = double(max(MSA_f_unique(:)));
weights = ones(B_f,1);  % re-weighting filtered MSA with threshold x = 1
MSA_id = sprintf('%s-N_%g-B_%g-x_1', dataID,N_f,B_f);


%% Correlation Compression
[MSA_cc, idx_cc] = CC_MSA(MSA_f_unique, MSA_id, N_f, B_f, q, weights, num_MI, ...
  numWorker, outputPath, NoLoad);
N_cc = numel(idx_cc);
B_cc = B_f;
CC_id = sprintf('CC-MI_%g-N_%g',num_MI,N_cc);
filename_CC_full = fullfile(outputPath, sprintf('%s--%s-MSA.mat',MSA_id,CC_id));
save(filename_CC_full, 'MSA_cc', 'idx_cc')


%% PLM
S = uint8(MSA_cc.'-1);


table_i_j_score = PLM_DCA(S,N_cc,B_cc,q,weights,lambda,numWorker);

% position in MSA_cc -> pos in MSA_f -> pos in original MSA (1-based)
table_i_j_score(1:2,:) = idx_f(idx_cc(table_i_j_score(1:2,:)));


%% save to file
DCA_id = sprintf('%s--PLM-l_%g',CC_id,lambda);

filename_score = sprintf('%s--%s.mat', MSA_id, DCA_id);
filename_score_full = fullfile(outputPath, filename_score);

fprintf('Saving to file ...\n')
tic
save(filename_score_full, 'table_i_j_score', '-v7.3');
time_save = toc;
fprintf('\tFinished in %.2f s.\n', time_save);


end
