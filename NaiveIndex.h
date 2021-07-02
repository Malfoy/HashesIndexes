#ifndef Index
#define Index

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>


typedef int8_t int8;
typedef uint8_t uint8;
typedef int32_t int32;
typedef uint32_t uint32;
typedef int64_t int64;
typedef uint64_t uint64;



using namespace std;



/*Structure of score which contains a trigger number of hits between query and
indexed genome and will be injected in Scores vector*/
struct score_strct {
uint32 genomeNumber;
double jaccardIndex;
};



class NaiveIndex {
public:
//attributes
int64 nb_minimizer;
int64 nb_genomes;
int64 bit_to_keep_minimizer;
int64 kmerSize;
vector<vector<uint8> > matrix;
/*Column=Partitions=Buckets line=genomesint64_t number_indexes;*/

//constructor
NaiveIndex(int64 bucketsNumber);

//methods
void add_sequence(const string& theSequence);
vector<score_strct> query_sequence(const string& sequenceSearched, int acceptanceTreshold); // the integral part represent the genome and decimals represent the Jaccard index
void add_sketch(const vector<uint8>& sketchToAdd);
vector<uint8> compute_sketch(const string& sequenceStr,const int kmerSize);
uint8 kmer_hasher(const string& str);
uint64 get_bucket(int64 primaryHash);
uint8 get_hash(int64 primaryHash);
};


#endif
