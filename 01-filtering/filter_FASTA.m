% Copyright (c) 2017 Chen-Yi Gao
%
% LICENSE
% ===
% See 'LICENSE.txt' in the outermost folder
%
% Description
% ===
% This function implements filtering based on FASTA file. See `filter_MSA` for
% details.
%
% HISTORY
% ===
% - 2017-10-24  v1.1
%   - add check on filename
%
% - 2017-10-17  v1
%   - adapted from `01filter.cpp`

function [MSA_f, idx_f] = filter_FASTA(filename, letter_N_max, MAF_min)

checkFilename(filename)

% search path
if exist('fasta2matrix_mex','file') ~= 3
  addpath(genpath(pwd))
end

fprintf('Reading FASTA file ...\n');
tic
% [MSA, len_seq, num_seq, num_dat] = fasta2matrix_mex(filename);
[MSA, len_seq, num_seq, ~] = fasta2matrix_mex(filename);
time_readFASTA = toc;
fprintf('\tFinished in %.2f s.\n', time_readFASTA);

MSA = MSA.';
q = double(max(MSA(:)));
B = double(num_seq);
N = double(len_seq);

[MSA_f, idx_f] = filter_MSA(MSA, B, N, q, letter_N_max, MAF_min);

end










% check if filename can be handled by `fasta2matrix_mex`
function checkFilename(filename)

if ~ischar(filename)
  error('filter_FASTA:filename', ...
    'filename should be provided as a char vector.\n  e.g. ''path/to/file''')
end

if contains(filename, '~')
  error('filter_FASTA:filename', ...
    '`~` is not supported, use `getenv(''HOME'')` instead.')
end

[status,values] = fileattrib(filename);
if ~status
  error('filter_FASTA:filename', ...
    '`%s` does not exist.', filename)
end
if ~(values.UserRead)
  error('filter_FASTA:filename', ...
    'Permission denied: `%s` can not be read.', filename)
end

end
