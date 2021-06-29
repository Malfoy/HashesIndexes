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



hash<std::string> kmerHasher;



NaiveIndex::NaiveIndex(int64_t n) //index constructor
{
  matrix.resize(n);
  nb_minimizer=n;
  bit_to_keep_minimizer=log2(n);
  nb_genomes=0;
  kmerSize=31;
};



uint64_t NaiveIndex::get_bucket(int64_t h)
{
  return h>>(64-bit_to_keep_minimizer); //bit shift to know bucket, first bits, most significant bit (MSB)
}



uint8_t NaiveIndex::get_hash(int64_t h)
{
  return h%256;  //to obtain the true hash, last bits, least significant bit (LSB)
}



vector<uint8_t> NaiveIndex::compute_sketch(const string& sequenceStr,const int kmerSize)
{
  vector<uint8_t> result(nb_minimizer,255); // the resulting vector

  // a For loop to hash every kmer of a sequence, distribute them in different
  // bucket with their MSB on one byte (255 buckets) and record some LSB number
  // of the minimal hash least

  int sequenceSize(sequenceStr.size()), counter(0);

  for (int position = 0 ; position < (sequenceSize - kmerSize) ; position++)
  {
    string kmerToHash(sequenceStr.substr (position,kmerSize));
    int64_t hashOfKmer(kmerHasher(kmerToHash));
    if(result[get_bucket(hashOfKmer)] > get_hash(hashOfKmer))
    {
      result[get_bucket(hashOfKmer)] =  get_hash(hashOfKmer);
    }
  }

  return result;
}
