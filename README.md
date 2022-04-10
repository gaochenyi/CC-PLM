# CC-PLM #

In each directory, there is a `README.md` concerning files therein. Please check them for details.

CC-PLM is a pilot version of correlation-compressed direct coupling analysis (CC-DCA) with the DCA flavour being PLM. All its components were revised to be user-friendly, however, there may be some ambiguities. Please let me know if any exists.

CC-PLM contains an optimized implementation of PLM for DCA (in directory `PLM-DCA`).

CC-PLM/PLM-DCA was developed as a part of an academic project.
**If you use CC-PLM/PLM-DCA (whether in whole or in part and whether modified or as is) for your own research, cite the following paper:**

- [Chen-Yi Gao, Hai-Jun Zhou and Erik Aurell, Correlation-Compressed Direct-Coupling Analysis, Phys. Rev. E 98, 032407.][link]

[link]: https://doi.org/10.1103/PhysRevE.98.032407

# How to reproduce our results? #

First, run `mexAll_CC_PLM` to compile all the MEX files. Then use the following snippets.

## CC-PLM ##

```matlab
fastafile = [getenv('HOME') '/loci/data/maela3K.fasta'];
dataID = 'maela3K';
outputPath = [getenv('HOME') '/loci/data'];
NoLoad = false;
numWorker = 56;

num_MI = 3e4;
lambda = 0.1;

paper_CC_PLM_DCA(fastafile,dataID, outputPath,NoLoad, numWorker, ...
    num_MI,lambda)
```

## PLM ##

``` matlab
fastafile = [getenv('HOME') '/loci/data/maela3K.fasta'];
dataID = 'maela3K';
outputPath = [getenv('HOME') '/loci/data'];
NoLoad = false;
numWorker = 56;

lambda = 0.1;
optTolOld = 1e-3;
optTolNew = 1e-3;
PLM_out_path = [getenv('HOME') '/loci/data/Potts'];

paper_PLM_DCA(fastafile,dataID, outputPath,NoLoad, numWorker, ...
    lambda,optTolOld,optTolNew,PLM_out_path)
```

## A summary about input arguments ##

These two `paper_*` functions are just for demonstration and reproduction of our results. For more detailed descriptions, please check `README.md` in subdirectories and comments along with the code in files.

- `NoLoad` and `outputPath`
  - At present, we have neither empirical nor theoretic rule to determine the hyper-parameter for PLM (`lambda`) and the hyper-parameter for CC (`num_MI`). Thus they must be specified by the user. One may want to run the code with different setting of these parameters, then part of the procedure shared by all runs can be stored in the first time and loaded later to avoid repeated calculation or processing of file. This feature is enabled by setting `NoLoad = false`, while `NoLoad = true` means always starting from the original FASTA file.
  - `outputPath` specifies the path for the output and intermediate results. Filename for output will be reported at the end.

| name      | meaning                                      |
| --------- | -------------------------------------------- |
| fastafile | full path for the FASTA file (See `README.md` in `01-filtering` for the limitation of the FASTA-handling function) |
| dataID    | identifier for data, which will be used for the filename of output |
| numWorker | number of workers for parallel computation   |
| lambda    | strength of the $l_2$ regularization for PLM |

### CC-PLM ###

`num_MI` specifies how many top correlations are used in the correlation-guided compression (CC) procedure.

### PLM ###

`paper_PLM_DCA` uses a modified version of `PLM_DCA`, `PLM_DCA_file`, which is dedicated to very large systems and allows resumption from previous partial run and making parameters more optimal.

`optTolNew` is the parameter for optimality of optimization involved in PLM.

By setting `optTolOld` equal to `optTolNew`, one can resume previous partial run.
For the first run, previous partial run is empty.

By setting `optTolOld = 1e-3` and `optTolNew = 1e-5`, one continues to optimize Potts parameters---this program tries to fetch Potts parameters from files, if any, associated with `1e-3`, and then continues the optimization with the stricter condition of optimality associated with `1e-5`.

FYI, it takes 14 days with `1e-3` on a 56-core server for an 81506-loci system.
It is estimated to take 10 more days to reach the stricter condition associated with `1e-5`.

Check the comment embedded in `PLM_DCA_file.m` for details.
