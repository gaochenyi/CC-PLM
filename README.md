# CC-PLM #

CC-PLM is a pilot version of correlation-compressed direct coupling analysis (CC-DCA) with the DCA flavour being PLM. All its components were revised to be user-friendly, however, there may be some ambiguities. Please let me know if any.

CC-PLM contains an optimized implementation of PLM for DCA (in directory `PLM-DCA`).

CC-PLM was developed as a part of an academic project. If possible, please cite the following paper:

  - [Chen-Yi Gao, Hai-Jun Zhou and Erik Aurell. Correlation-compressed Direct Coupling Analysis. arXiv:1710.04819][our paper]

[our paper]: https://arxiv.org/abs/1710.04819


# How to reproduce our results? #

First, run `mexAll_CC_PLM` to compile all the MEX files. Then use the following snippets.

## CC-PLM ##

```matlab
fastafile = [getenv('HOME') '/loci/data/maela3K.fasta'];
dataID = 'maela3K';
lambda = 0.1;
outputPath = [getenv('HOME') '/loci/data'];
NoLoad = false;
numWorker = 56;

num_MI = 3e4;
paper_CC_PLM_DCA(fastafile, dataID, num_MI, lambda, numWorker, outputPath, NoLoad)
```

## PLM ##

``` matlab
fastafile = [getenv('HOME') '/loci/data/maela3K.fasta'];
dataID = 'maela3K';
lambda = 0.1;
outputPath = [getenv('HOME') '/loci/data'];
NoLoad = false;
numWorker = 56;

PLM_out_path = [getenv('HOME') '/loci/data/Potts'];
optTolOld = 1e-3;
optTolNew = 1e-3;
paper_PLM_DCA(fastafile,dataID,lambda,outputPath,NoLoad,optTolOld,optTolNew,PLM_out_path,numWorker)
```
