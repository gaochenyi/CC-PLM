# Summary

This directory contains programs necessary for filtering:

- `filter_MSA` performs filtering based on MSA matrix, rows of which are sequences.
- `filter_FASTA` constructs a MSA matrix from a [FASTA][] file, then use `filter_MSA` performs the filtering procedure.
- `mex_fasta` compiles the MEX file for reading FASTA file.
- `test.fasta` is an example FASTA file.
- The directory `function` contains supporting functions.

[FASTA]: https://www.ncbi.nlm.nih.gov/BLAST/fasta.shtml

# FASTA format

At present only a small subset of [FASTA format][FASTA] is supported. The specification is:

1. For each sequence, a single-line comment precedes lines of data. The comment line starts with `>`. Consecutive data lines will be concatenated into one long data line.
2. Lower-case letters are accepted and are mapped into upper-case.
3. The nucleic acid codes supported are: NACGT
