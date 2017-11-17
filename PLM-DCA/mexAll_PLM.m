function mexAll_PLM

here = pwd;

fprintf(['\n' ...
  'PLM for Potts model\n' ...
  '===================\n'])

if exist('function/compiled', 'dir') ~= 7
  mkdir('function/compiled')
end
fprintf('Compiling `g_r_mex_v2.cpp` ...\n')
mex -silent -largeArrayDims CXXFLAGS='$CXXFLAGS -std=c++11' ...
  -outdir function/compiled function/mex/g_r_mex_v2.cpp


% minFunc (third party)
fprintf(['\n' ...
  'Third Party Software\n' ...
  '--------------------\n'])
cd ./thirdparty
% MEX file in minFunc uses MATLAB Version 7.2 array-handling API
if exist('minFunc/compiled', 'dir') ~= 7
  mkdir('minFunc/compiled')
end
fprintf('Compiling minFunc files...\n');
mex -compatibleArrayDims -outdir minFunc/compiled minFunc/mex/mcholC.c
mex -compatibleArrayDims -outdir minFunc/compiled minFunc/mex/lbfgsC.c
mex -compatibleArrayDims -outdir minFunc/compiled minFunc/mex/lbfgsAddC.c
mex -compatibleArrayDims -outdir minFunc/compiled minFunc/mex/lbfgsProdC.c
cd(here)

end
