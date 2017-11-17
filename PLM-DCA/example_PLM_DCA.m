clear

N = 20;
B = 3000;
q = 4;  % nucleic acid case
MSA = randi(q,N,B);

%% PLM
S = uint8(MSA-1);   % [1,q] -> [0,q-1]
weights = ones(B,1);
lambda = 0.1;
numWorker = 2;

table_i_j_score = PLM_DCA(S,N,B,q,weights,lambda,numWorker);
