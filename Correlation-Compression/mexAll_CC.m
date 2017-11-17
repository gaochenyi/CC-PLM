fprintf(['\n' ...
  'Correlation Compression\n' ...
  '=======================\n'])

if exist('function/compiled', 'dir') ~= 7
  mkdir('function/compiled')
end

fprintf('Compiling `calc_f2_w_mex_uint8.cpp` ...\n')
mex -silent -largeArrayDims CXXFLAGS='$CXXFLAGS -std=c++11' ...
  -outdir function/compiled function/mex/calc_f2_w_mex_uint8.cpp
