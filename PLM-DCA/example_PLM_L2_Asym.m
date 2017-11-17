clear

N = 31;
B = 3001;
q = 11;
MSA = randi(q,N,B);

%% PLM
S = uint8(MSA-1);
weights = ones(B,1);
lambda = 1;
lambdas = [lambda lambda/2];  % to be compatible with symmetric version

numWorker = 2;

skip = false; % Built-in checks of g_r_mex_v2 only cause small overhead.


%% see `minFunc.m` for details
options.Display = 'off';    % 'off' to not display progress information
options.progTol = -1;       % stop only when 1st-order optimality reaches (controlled by `optTol`)
options.optTol  = 1e-5;     % smaller optTol -> closer to optimal point, but more iterations
options.useMEX  = true;     % Tests show that MEX boosts by 15%
options.Method  = 'lbfgs';  % L-BFGS (2nd-order method)
% number of corrections to store in memory, used to construct a approximation of
% Hessian, more corrections result in faster convergence but use more memory
options.Corr    = 100;      % (default: 100)

h_and_J = PLM_L2_Asym(S,N,B,q,weights,lambdas,skip,options,numWorker);
