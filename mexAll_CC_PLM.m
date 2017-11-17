function mexAll_CC_PLM

here = pwd;

cd 01-filtering
mex_fasta
cd(here)

cd Correlation-Compression
mexAll_CC
cd(here)

cd PLM-DCA
mexAll_PLM
cd(here)

end
