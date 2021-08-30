#ifndef INDEX
#define INDEX
#define NDEBUG  //ASSERT ON TESTTABLE.CPP

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <bit>
#include <bitset>
#include <stdint.h>
#include <utility>
#include <tr1/unordered_map>

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





/* Produce an index matrix containing a sketch by genome line with the function
   to add sketch and to query a sequence among the genome with a score based
   on jaccard index. */
class NaiveIndex {
public:


//~~Attributes~~
vector<int> goodGenomeVector;
uint64 nb_minimizer;
uint64 nb_genomes;
uint64 decimal_lsb;
uint64 number_of_LSB;
uint64 bit_to_keep_minimizer;
uint64 kmerSize;
vector<vector<uint8> > matrix; //Column==Buckets line==genomes


//~~Constructor~~
NaiveIndex(uint64 bucketsNumber,uint16 decimal_for_lsb = 256,uint16 nbgenomes = 0);


//~~Methods~~
//    ~~public~~
void index_sequences_from_fasta(const string& fileName);
vector<vector<long double>> get_all_queries_scores_for_index(vector<vector<long double>> ultimVector, const string& fileName);
vector<int> get_good_genome_vector();
vector<long double> query_sequence(string sequenceSearched);
vector<pair<long double,uint16>> sort_scores(vector<long double> allScoresVector, long double thresholdJaccard = 0);
void show_sorted_scores(vector<pair<long double,uint16>> sortedScoresVector, uint howManyScoresToShow = 0);

//    ~~private~~
string get_line_fasta_for_naive(ifstream* partToExamine);
vector<uint8> compute_sketch(string sequenceStr);
void add_sketch(const vector<uint8>& sketchToAdd);

  string get_complement_or_not(string sequenceToComplement);
  uint64 get_bucket(uint64 primaryHash);
  uint8 get_hash(uint64 primaryHash);

    uint8 get_minhash(uint64 primaryHash);

//    ~~alternative~~
uint8 get_hyperloglog(uint64 primaryHash);
uint64 xs(uint64 y);
uint8 get_popcount(uint64 primaryHash);
uint8 get_hyper_minhashX(uint64 primaryHash, uint8 hyperMinhashNumber = 62);
uint8 get_double_hyperloglog(uint64 primaryHash);
};

#endif
