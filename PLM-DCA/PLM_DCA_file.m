% Copyright (c) 2017 Chen-Yi Gao
% 
% LICENSE
% ===
% The MIT License
% 
% 
% DESCRIPTION
% ===
% Designed to be used on very large systems (e.g., 80000-node system).
% 
% This function try to boost the optimization (with optimal tolerance
% `optTolNew`) by results obtained with optimal tolerance `optTolOld` (if they
% exist). Thus,
% 
% 1. when optTolOld = optTolNew, previous partial run will be resumed;
% 2. when optTolOld > optTolNew), parameters will be more optimal.
% 
% For very large systems, DCA will take a few days with default optimal
% tolerance (1e-5). Relaxing the optimal tolerance will lead to shorter runtime.
% But maybe later one want to see whether results change when Potts parameters
% are more optimal---this is the scenario 2.
% 
% Since DCA for large system will take a long time (even with relaxed optimal
% tolerance 1e-3), it may quits accidentally, due to HPC issue or something
% else. One then want to resume the previous partial run---this is scenario 1.
% 
% Initial run means that previous partial run is empty---this is also
% scenario 1.
% 
% My test suggest that it is acceptable that `optTol=1e-3` but not further
% relaxed. (For small system, default value of `optTol` is 1e-5.)
% 
% 
% HISTORY
% ===
% - 2017-11-17  v1.1
%   - PLM_L2_Asym_resume -> PLM_L2_Asym_file
% 
% - 2017-10-31  v1

function table_i_j_score = PLM_DCA_file(S,N,B,q,weights,lambda,numWorker, ...
  optTolOld,optTolNew,PLM_out_path,filePrefix)

% search path
addpath(genpath(pwd))


%% PLM
% see `minFunc.m` for details
options.Display = 'off';    % 'off' to not display progress information
options.progTol = -0;       % stop only when 1st-order optimality reaches (controlled by `optTol`)
options.optTol  = optTolNew;% smaller optTol -> closer to optimal point, but more iterations
options.useMEX  = true;     % Tests show that MEX boosts by 15% for N~8e4, B~3e3
options.Method  = 'lbfgs';  % L-BFGS (2nd-order method)
% number of corrections to store in memory, used to construct a approximation of
% Hessian, more corrections result in faster convergence but use more memory
options.Corr    = 100;      % (default: 100)

lambdas = [lambda lambda/2];  % Every J_{ij}(a,b) counts twice in the asymmetric version.
skip = false;
LoadIP = true;
SaveFP = true;
filePrefixLoad = sprintf('%s-opt_%g',filePrefix,optTolOld);
filePrefixSave = sprintf('%s-opt_%g',filePrefix,optTolNew);
h_and_J = PLM_L2_Asym_file(S,N,B,q,weights,lambdas,skip,options,numWorker, ...
  LoadIP,SaveFP,PLM_out_path,filePrefixLoad,filePrefixSave);


%% Scoring couplings
fprintf('Scoring the coupling ...\n')
timer = tic;
table_i_j_score = score_coupling_L2_no_gap(h_and_J,q,N);
time = toc(timer);
fprintf('\tFinished in %.2f s.\n', time);


end
