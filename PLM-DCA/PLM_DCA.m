% Copyright (c) 2017 Chen-Yi Gao
% 
% LICENSE
% ===
% The MIT License
% 
% HISTORY
% ===
% - 2017-10-24  v1

function table_i_j_score = PLM_DCA(S,N,B,q,weights,lambda,numWorker)

% search path
addpath(genpath(pwd))


%% PLM
% see `minFunc.m` for details
options.Display = 'off';    % 'off' to not display progress information
options.progTol = -0;       % stop only when 1st-order optimality reaches (controlled by `optTol`)
options.optTol  = 1e-5;     % smaller optTol -> closer to optimal point, but more iterations
options.useMEX  = true;     % Tests show that MEX boosts by 15% for N~8e4, B~3e3
options.Method  = 'lbfgs';  % L-BFGS (2nd-order method)
% number of corrections to store in memory, used to construct a approximation of
% Hessian, more corrections result in faster convergence but use more memory
options.Corr    = 100;      % (default: 100)

lambdas = [lambda lambda/2];  % Every J_{ij}(a,b) counts twice in the asymmetric version.
skip = false;
h_and_J = PLM_L2_Asym(S,N,B,q,weights,lambdas,skip,options,numWorker);


%% Scoring couplings
fprintf('Scoring the coupling ...\n')
timer = tic;
table_i_j_score = score_coupling_L2_no_gap(h_and_J,q,N);
time = toc(timer);
fprintf('\tFinished in %.2f s.\n', time);


end
