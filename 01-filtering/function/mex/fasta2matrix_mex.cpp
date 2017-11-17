/**
 * Copyright (c) 2017 Chen-Yi Gao
 *
 * LICENSE
 * ===
 * See 'LICENSE.txt' in the outermost folder
 *
 *
 * DESCRIPTION
 * ===
 * MEX wrapper of `fasta2matrix.cpp`
 *
 * This function reads FASTA file for DNA sequences. Supported nucleic acid
 * codes are: N and ATGC. Assume that size_t is equal to uint64_t
 *
 *
 * HISTORY
 * ===
 * - 2017-10-17  v2
 *   - NACGT-01234 --> NACGT-12345
 *
 * - 2017-09-25  v1
 *  - initial draft
 */

#include <iostream> // test
#include <fstream>
#include <string>
#include <vector>
// #include <cstdint>
#include "mex.h"

using namespace std;


// typedef int8_T Type;
// #define Type_CLASS_ID mxINT8_CLASS
typedef uint8_T Type;                   // type of number representation
#define Type_CLASS_ID mxUINT8_CLASS



// letter to number: NACGT -> 12345
template <class T>
inline T let2num(const char& let) {
  T num;
  switch(let) {
    case 'N' :
    case 'n' : num = 1; break;
    case 'A' :
    case 'a' : num = 2; break;
    case 'C' :
    case 'c' : num = 3; break;
    case 'G' :
    case 'g' : num = 4; break;
    case 'T' :
    case 't' : num = 5; break;
    default :
      // error
      mexErrMsgIdAndTxt("fasta2matrix:let2num",
        "Unsupported letter: %c\n", let);
  }
  return num;
}



// string sequence to number sequence
inline vector<Type> str2vec(const string& seq_str) {
  vector<Type> seq_num;
  seq_num.reserve(seq_str.length());
  for (const auto &c : seq_str) {
    seq_num.emplace_back(let2num<Type>(c));
  }
  return seq_num;
}



// given filename of MSA in FASTA, return the corresponding numeric matrix
inline void fasta2mxArray(
  const char* filename_in,
  mxArray* &pm_MSA,
  mxArray* &pm_len_seq, mxArray* &pm_num_seq, mxArray* &pm_num_dat)
{
  ifstream fin(filename_in);
  if(!fin) {
    mexErrMsgIdAndTxt("fasta2matrix:file",
      "Could not read file '%s'.\n", filename_in);
  }

  /* FASTA -> number sequences */
  size_t num_seq = 0; // number of sequences processed
  size_t num_dat = 0; // number of data lines processed
  typedef vector<Type> seqType;
  vector<seqType> msa_num;
  for(string line_str; getline(fin, line_str); ) {
    // skip empty lines
    if (line_str.length() == 0) {
      continue;
    }

    // a comment line precedes a sequence
    if (line_str[0] == '>') {
      num_seq++;
      msa_num.emplace_back(seqType());
      continue;
    }

    if (num_seq == 0) {
      // first data line without comment header
      mexErrMsgIdAndTxt("fasta2matrix:FASTA",
        "FASTA file is illegal---no comment precedes the first data line.\n");
    }
    else {
      // merge consecutive data lines to one sequence
      num_dat++;
      auto &seq_num = msa_num[num_seq-1];
      const auto &line_num = str2vec(line_str);
      seq_num.reserve(seq_num.size() + line_num.size());
      seq_num.insert(seq_num.end(), line_num.begin(), line_num.end());
    }
  }
  // check EOF
  if (!fin.eof()) {
    mexErrMsgIdAndTxt("fasta2matrix:file",
      "File parsing stops before reaching EOF.\n");
  }

  /* check */
  if (num_seq == 0) {
    mexErrMsgIdAndTxt("fasta2matrix:FASTA:NoSequence",
      "'%s' contains no sequence.\n", filename_in);
  }

  const size_t len_seq = msa_num[0].size(); // msa_num.size() >= 1

  if (len_seq == 0) {
    mexErrMsgIdAndTxt("fasta2matrix:FASTA:FirstSequenceVoid",
      "The length of the first sequence is 0.");
  }

  // check if MSA is legal
  for (size_t i = 1; i < num_seq; i++) {
    if (msa_num[i].size() != len_seq) {
      mexErrMsgIdAndTxt("fasta2matrix:MSA",
        "The length of sequence %ld doesn't match that of sequence 1.", i+1);
    }
  }

  // // test
  // cout << num_seq << " sequences\n";
  // cout << num_dat << " data lines\n";
  // cout << "length = " << len_seq << '\n';
  // for (const auto &i : msa_num) {
  //   for (const auto &j : i) {
  //     cout << int(j);
  //   }
  //   cout << '\n';
  // }

  // copy MSA to mxArray
  pm_MSA = mxCreateNumericMatrix(len_seq, num_seq, Type_CLASS_ID, mxREAL);
  Type* pr = (Type*) mxGetData(pm_MSA);
  size_t idx = 0;
  for (const auto &seq : msa_num) {
    for (const auto &let : seq) {
      pr[idx] = let;
      idx++;
    }
  }

  // set len_seq, num_seq, num_dat
  pm_len_seq = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
  *((size_t*) mxGetData(pm_len_seq)) = len_seq;

  pm_num_seq = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
  *((size_t*) mxGetData(pm_num_seq)) = num_seq;

  pm_num_dat = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
  *((size_t*) mxGetData(pm_num_dat)) = num_dat;
}


/* [MSA, len_seq, num_seq, num_dat] = fasta2matrix_mex(filename) */
void mexFunction(
  int nlhs, mxArray *plhs[],
  int nrhs, const mxArray *prhs[])
{
  // validation
  if (nrhs != 1) {
    mexErrMsgTxt("This function accepts 1 argument.");
  }
  if (nlhs > 4) {
    mexErrMsgTxt("This function produces at most 4 outputs.");
  }
  const mxArray* pfilename = prhs[0];
  if (!mxIsChar(pfilename)) {
    mexErrMsgTxt("`filename` should be provided as a string.");
  }

  char* pc = mxArrayToString(pfilename);
  if (pc == NULL) {
    mexErrMsgTxt("`mxArrayToString` failed.");
  }

  // real work
  fasta2mxArray(pc, plhs[0], plhs[1], plhs[2], plhs[3]);

  mxFree(pc);
}