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
  vector<vector<uint8_t>> matrix;//Column=Partitions line=genomes

  NaiveIndex NaiveIndex(int64_t n);
  void add_sequence(const string& str);
  vector<double> query_sequence(const string& str);

  void add_sketch(const vector<uint8_t>& V);
  vector<uint8_t> compute_sketch(const string& str);
  uint8_t hash(const string& str);

};


#endif
