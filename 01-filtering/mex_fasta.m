fprintf(['\n' ...
  'Filtering procedure\n' ...
  '===================\n'])

if exist('function/compiled', 'dir') ~= 7
  mkdir('function/compiled')
end

fprintf('Compiling `fasta2matrix_mex.cpp` ...\n')
mex -silent -largeArrayDims CXXFLAGS='$CXXFLAGS -std=c++11' ...
  -outdir function/compiled function/mex/fasta2matrix_mex.cpp
