% Copyright (c) 2017 Chen-Yi Gao
%
% LICENSE
% ===
% See 'LICENSE.txt' in the outermost folder
%
% DESCRIPTION
% ===
% Implementation of all operations needed on each locus for filtering.
%
% Assumption
% ===
% `data`, a column vector of type uint8, contains samples of a locus. Supported
% letters are: NACGT, which are represented as 12345.
%
% **No sophiscated check on input data!**
%
% memo for author: See p88 of notebook for more comprehensive procedure.
%
% HISTORY
% ===
% - 2017-10-17  v1
%   - adapted from `01filter.cpp`

function [data_new_rep, label] = filter_locus(...
  data, B, q, ...
  N_counter_max, MAF_min)


%% simple check
if q ~= 5 || size(data,1) ~= numel(data) || ~isa(data, 'uint8')
  error('Assumption fails.')
end

%% real work
counter = calc_counter(data, B, q);

% gap-rich
gapRich = (counter(1) > N_counter_max);

% Ignoring N, find 1st/2nd/3rd/4th most common letter. in sidx: 1234 --> ACGT
[counter_sort,sidx] = sort(counter(2:5),'descend');

newRep = uint8((1:5).');
newRep(sidx+1) = uint8(2:5);

% multi-allelic
multiAllelic = (counter_sort(3) > 0);

% MAF-low: MAF < MAF_min; MAF = minor / (minor + major)
minorPoor = (counter_sort(2) / (counter_sort(2) + counter_sort(1)) < MAF_min);

% bit-like label, 0~7
label = gapRich*4 + minorPoor*2 + multiAllelic;

% N, 1st, 2nd, 3rd, 4th representation of data (as 12345)
data_new_rep = newRep(data);

end
