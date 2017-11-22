% The customized version to produce results in our paper

% function paper_PLM_DCA(fastafile,dataID,lambda, outputPath,NoLoad, optTolOld,optTolNew,PLM_out_path,numWorker)
function paper_PLM_DCA(fastafile,dataID, outputPath,NoLoad, numWorker, ...
  lambda, optTolOld,optTolNew,PLM_out_path)

%% check with little overhead

% dataID
if ~ischar(dataID)
  error('`dataID` should be provided as a char vector.')
end

% lambda
if lambda < 0.0
  error('lambda should be non-negative.')
end

% check paths
checkPath(outputPath)
checkPath(PLM_out_path)

% NoLoad
if ~islogical(NoLoad)
  error('`NoLoad` should be provided as logical.')
end

% optTol
if optTolNew <= 0
  error('`optTol` need to be positive.')
elseif optTolNew > 1e-3
  error('It is not recommended to set optTol larger than 1e-3.')
end


%% some preparation

% If no pool exists, a new one is created.
poolobj = gcp('nocreate');
if isempty(poolobj)
  parpool(numWorker);
end

% search path
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


%% PLM
DCA_id = sprintf('PLM-l_%g-opt_%g', lambda, optTolNew);
filePrefix = sprintf('%s--PLM-parameters-l_%g',MSA_id,lambda);

S_f = uint8(MSA_f_unique.'-1);
table_i_j_score = PLM_DCA_file(S_f,N_f,B_f,q,weights,lambda,numWorker, ...
  optTolOld,optTolNew,PLM_out_path,filePrefix);
% pos in MSA_f -> pos in original MSA (1-based)
table_i_j_score(1:2,:) = idx_f(table_i_j_score(1:2,:));


%% save to file

filename_score = sprintf('%s--%s.mat', MSA_id, DCA_id);
filename_score_full = fullfile(outputPath, filename_score);

fprintf('Saving to file ...\n')
tic
save(filename_score_full, 'table_i_j_score', '-v7.3');
time_save = toc;
fprintf('\tFinished in %.2f s.\n', time_save);

fprintf('Full path to the output file is \n\n\t%s\n\n',filename_score_full)

end


% check whether path is legal
function checkPath(Path)

if ~ischar(Path)
  error('Path should be provided as a char vector.')
end
if exist(Path,'dir') ~= 7
  error('The folder `%s` does not exist.', Path)
end

end
