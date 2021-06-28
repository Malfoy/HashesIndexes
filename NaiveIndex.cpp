#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <functional>
#include <mutex>
#include <string>
#include <stdint.h>
#include <math.h>
#include "NaiveIndex.h"



using namespace std;



hash<std::string> str_hash;



NaiveIndex::NaiveIndex(int64_t n){
  matrix.resize(n);
  nb_minimizer=n;
  bit_to_keep_minimizer=log2(n);
  nb_genomes=0;
};



uint64_t NaiveIndex::get_bucket(int64_t h){
  return h>>(64-bit_to_keep_minimizer);
}



uint8_t NaiveIndex::get_hash(int64_t h){
  return h%256;
}



vector<uint8_t> NaiveIndex::compute_sketch(const string& str){
  vector<uint8_t> result(nb_minimizer,255);


  return result;
}
