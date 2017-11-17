% Copyright (c) 2017 Chen-Yi Gao
% 
% LICENSE
% ===
% The MIT License
% 
% DESCRIPTION
% ===
% Given B weighted samples, this function performs the asymmetric version of
% L2-regularized pseudo-likelihood maximization (PLM) for q-state Potts model.
% The inferred parameters are transformed to Ising gauge (zero-sum gauge).
%
% INPUT
% ===
% S  (uint8)  [0,q-1], columns as sequences/samples/configurations; gap state should be mapped to 0
% N           length of sequences (number of nodes/spins)
% B           number of sequences/samples/configurations
% q           number of possible states
% weights     [0,1]  sequences can be weigted
% lambdas     [lambda_h lambda_J]
% skip        non-zero to skip built-in check of `g_r_mex_v2`
% options     passed to `minFunc`
% numWorker   number of workers used in parfor
%
% OUTPUT
% ===
% h_and_J(:,r) represents [h_r(:); J_r(:)], where h_r and J_r are in Ising gauge
%
% HISTORY
% ===
% - 2017-11-16  v2.1
%   - add check for limitation
% 
% - 2017-10-26  v2
%   - PLM and gauge transformation are now seperated (to boost parfor).
% 
% - 2017-10-23  v1.1
%   - Temporary variable `wr` is eliminated. (small speed-up)
% 
% - 2017-10-19  v1.0
%   - initial draft

function h_and_J = PLM_L2_Asym(S,N,B,q,weights,lambdas,skip,options,numWorker)

if nargin ~= 9
  error('Not enough input arguments.')
end

%% check with very little overhead
% type
if ~isa(S,'uint8') ...
    || ~isa(N, 'double') || ~isa(B, 'double') || ~isa(q, 'double') ...
    || ~isa(weights, 'double') || ~isa(lambdas, 'double') ...
    || ~isa(numWorker, 'double')
  error(...
    ['%s:\n' ...
    '   uint8:    S\n' ...
    '  double:    N, B, q, r, weights, lambdas, numWorker'], ...
    'Requirement on type');
end

% orientation of S
if size(S,1) ~= N || size(S,2) ~= B
  error('`S` should be a N-by-B matrix.')
end

% dimension of weights
if numel(weights) ~= B
  error('weights should contains B numbers.')
end

% check whether N/B/q is a integer
if round(N) ~= N || round(B) ~= B || round(q) ~= q ...
    || round(numWorker) ~= numWorker
  error('N, B, q and numWorker should be integers.')
end

% check if limitation is reached
if q > 256
  error('At most 256 states are supported.')
end
if options.useMEX && strcmp(options.Method, 'lbfgs') ...
    && options.Corr*(q + q*q*(N-1)) >= 2^31
  error('Limitation is reached; extra work needed; see `README.md`.')
end


%% check with acceptable overhead

% range of S
if double(max(S(:))) > q-1
  error('q possible states should be mapped to integers in [0,q-1].')
end

% range of weights
for b = 1:B
  if weights(b) < 0 || weights(b) > 1
    error('`weights` may not exceeds [0,1].')
  end
end


%% search path

if exist('min_g_r', 'file') ~= 2
  addpath(genpath(pwd))
end


%% real work

% If no pool exists, a new one is created when necessary
if numWorker > 1
  poolobj = gcp('nocreate');
  if isempty(poolobj)
    parpool(numWorker);
  end
end


B_eff = sum(weights);
h_and_J = zeros(q + q*q*(N-1), N);


%% PLM
fprintf('Performing L2-regularized PLM (asymmetric version) ...\n')
timer = tic;

if numWorker > 1
  parfor (r = 1:N, numWorker)
    h_and_J(:,r) = min_g_r( ...
      S, uint64(N),uint64(B),uint64(q), ...
      weights,B_eff,uint64(r),lambdas,skip,options);
  end
else
  for r = 1:N
    h_and_J(:,r) = min_g_r( ...
      S, uint64(N),uint64(B),uint64(q), ...
      weights,B_eff,uint64(r),lambdas,skip,options);
  end
end

time = toc(timer);
fprintf('\tFinished in %.2f s.\n', time);


%% Gauge Transformation
fprintf('Shifting parameters of Potts model to Ising gauge ...\n');
timer = tic;

if numWorker > 1
  parfor (r = 1:N, numWorker)
    h_and_J(:,r) = gauge_shift_Ising(h_and_J(:,r), q, N);
  end
else
  for r = 1:N
    h_and_J(:,r) = gauge_shift_Ising(h_and_J(:,r), q, N);
  end
end

time = toc(timer);
fprintf('\tFinished in %.2f s.\n', time);


end
