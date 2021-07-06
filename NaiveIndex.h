#ifndef Index
#define Index

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <bit>
#include <bitset>
#include <stdint.h>

// removal '_t'
typedef int8_t int8;
typedef uint8_t uint8;
typedef int16_t int16;
typedef uint16_t uint16;
typedef int32_t int32;
typedef uint32_t uint32;
typedef int64_t int64;
typedef uint64_t uint64;



using namespace std;



/*Structure of score containing a trigger number of hits between query and
   indexed genome and will be injected in Scores vector*/
struct score_strct {
        uint32 genomeNumber;
        double jaccardIndex;
};


/* Produce an index matrix containing a sketch by genome line with the function
   to add sketch and to query a sequence among the genome with a score based
   on jaccard index. */
class NaiveIndex {
public:

//attributes
int64 nb_minimizer;
int64 nb_genomes;
uint16 decimal_lsb;
int64 bit_to_keep_minimizer;
int64 kmerSize;
vector<vector<uint8> > matrix; //Column==Buckets line==genomes


//constructor
NaiveIndex(int64 bucketsNumber,uint16 decimal_for_lsb = 256);

//methods
void add_sequence(const string& theSequence);
vector<score_strct> query_sequence(const string& sequenceSearched, int acceptanceTreshold);
void add_sketch(const vector<uint8>& sketchToAdd);
vector<uint8> compute_sketch(const string& sequenceStr,const int kmerSize);
uint8 kmer_hasher(const string& str);
uint64 get_bucket(int64 primaryHash);
uint8 get_hash(int64 primaryHash);
uint8 get_minhash(int64 primaryHash);
uint8 get_hyperloglog(int64 primaryHash);
uint8 get_popcount(int64 primaryHash);
uint8 get_hyper_minhash62(int64 primaryHash);
uint8 get_hyper_minhash53(int64 primaryHash);
uint8 get_hyper_minhash44(int64 primaryHash);
uint8 get_hyper_minhash35(int64 primaryHash);
uint8 get_hyper_minhash26(int64 primaryHash);
};

class TestTable {
  
};
#endif
