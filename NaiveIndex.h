#ifndef Index
#define Index



#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>



using namespace std;


class NaiveIndex{
public:
      int64_t nb_minimizer;
      int64_t nb_genomes;
      int64_t bit_to_keep_minimizer;
      int64_t kmerSize;
      vector<vector<uint8_t>> matrix;//Column=Partitions line=genomesint64_t number_indexes;
      NaiveIndex(int64_t n);
      void add_sequence(const string& str);
      vector<double> query_sequence(const string& str);

      void add_sketch(const vector<uint8_t>& V);
      vector<uint8_t> compute_sketch(const string& sequenceStr,const int kmerSize);
      uint8_t hash_str(const string& str);
      uint64_t get_bucket(int64_t h);
      uint8_t get_hash(int64_t h);
};


#endif
