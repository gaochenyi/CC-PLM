% Copyright (c) 2017 Chen-Yi Gao
% 
% LICENSE
% ===
% The MIT License
% 
% DESCRIPTION
% ===
% minimization of g_r, with support of loading and/or saving results. When
% loading and saving are disabled, this function is essentially same as
% `min_g_r`.
% 
% INPUT
% ===
% |    name    | description                                               |
% | ---------- | --------------------------------------------------------- |
% | S (uint8)  | [0,q-1], columns as sequences/samples/configurations      |
% | N (uint64) | length of sequences (number of nodes/spins)               |
% | B (uint64) | number of sequences/samples/configurations                |
% | q (uint64) | number of possible states                                 |
% | weights    | sequences can be weigted                                  |
% | B_eff      | B_eff = sum(weights)                                      |
% | r (uint64) | node index (1-based)                                      |
% | lambdas    | [lambda_h lambda_J]                                       |
% | skip       | non-zero to skip built-in check of `g_r_mex_v2`           |
% | options    | passed to `minFunc`                                       |
% | LoadIP     | true to load initial point                                |
% | SaveFP     | true to save final point                                  |
% | filePath   | path of files                                             |
% | filePrefix | the filename will be constructed by [filePath filePrefix] |
% 
% OUTPUT
% ===
% r_h_and_J = [h_r(:); J_r(:)], where h_r and J_r are in Ising gauge.
% 
% HISTORY
% ===
% - 2017-11-16  v1
%   - adapted from `min_g_r.m`

function r_h_and_J = min_g_r_file(S,N,B,q,weights,B_eff,r,lambdas,skip,options, ...
  LoadIP,SaveFP,filePath,filePrefixLoad,filePrefixSave)

% simple check
if ~ischar(filePath) || ~ischar(filePrefixLoad) || ~ischar(filePrefixSave)
  error('path and filename should be provided as char vectors.')
end

if exist(filePath,'dir') ~= 7
  error('`%s` does not exist.', filePath);
end


% Load the result if it exists, otherwise start from zeros
r_h_and_J = 0;
filenameLoad = sprintf('%s--r-%d',filePrefixLoad,r);
filenameLoadFull = fullfile(filePath, [filenameLoad '.mat']);
if LoadIP && exist(filenameLoadFull, 'file') == 2
  load(filenameLoadFull, 'r_h_and_J');        % load initial point
  if size(r_h_and_J,1) ~= (q+q*q*(N-1)) ...   % check dimension of `r_h_and_J`
      || numel(r_h_and_J) ~= (q+q*q*(N-1))
    error('Dimension of `r_h_and_J` should be [q+q*q*(N-1), 1].')
  end
else
  r_h_and_J = zeros(q + q*q*(N-1), 1);
end


% minimization of g_r
if skip
  funObj = @(wr) g_r_mex_v2(S,N,B,q,weights,B_eff,r,wr,lambdas,'SkipCheckFlag');
else
  funObj = @(wr) g_r_mex_v2(S,N,B,q,weights,B_eff,r,wr,lambdas);
end
r_h_and_J = minFunc(funObj,r_h_and_J,options);


% save the result to file
if SaveFP
  filenameSave = sprintf('%s--r-%d',filePrefixSave,r);
  filenameSaveFull = fullfile(filePath, [filenameSave '.mat']);
  if numel(r_h_and_J)*8 > 2^31
    error('Extra work needed; see the embedded comment.')
  end
  save(filenameSaveFull,'r_h_and_J','options','-v6');
  % '-v6' limits the maximum size of any variable to 2^31 bytes. For larger
  % variables, use HDF5 without compression (by `h5write` and friends).
end


end
