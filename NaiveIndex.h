#ifndef Index
#define Index
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
int64 nb_minimizer;
int64 nb_genomes;
uint16 decimal_lsb;
int64 bit_to_keep_minimizer;
int64 kmerSize;
vector<vector<uint8> > matrix; //Column==Buckets line==genomes

//~~Constructor~~
NaiveIndex(int64 bucketsNumber,uint16 decimal_for_lsb = 256);

//~~Methods~~
void add_sequence(const string& theSequence);
vector<double> query_sequence(const string& sequenceSearched, double acceptanceTreshold = 0);
void add_sketch(const vector<uint8>& sketchToAdd);
vector<uint8> compute_sketch(const string& sequenceStr,const int kmerSize);
uint8 kmer_hasher(const string& str);
uint64 get_bucket(int64 primaryHash);
uint8 get_hash(int64 primaryHash);
uint8 get_minhash(int64 primaryHash);
uint8 get_hyperloglog(int64 primaryHash);
uint8 get_popcount(int64 primaryHash);
uint8 get_hyper_minhashX(int64 primaryHash, uint8 hyperMinhashNumber = 62);
uint8 get_double_hyperloglog(int64 primaryHash);
};



class TestTable {
public:

//~~Attributes~~
int64 kmerSize;
uint32 nbGenomes;

//~~Constructor~~
TestTable(uint32 genomeQuantityForTest);

//~~Methods~~
void record_sequence(const string& sequenceStr, const uint32 Genome);
vector<pair <uint32, uint32> > query_belonging_genome(const string& sequenceStr);
};



// class about a structure to evaluate results
class ComparisonMatrix {
public:

//~~Attribute~~
uint32 matrixHeight;

//~~Constructor~~
ComparisonMatrix();

//~~Methods~~
void add_result_vector(vector<double> resultVector);
void create_comparison();
void print_vector(vector<double> vectorToPrint);
void show_the_matrix();
};



// class to extract sequences from fasta file into string type
class FastaParser{
public:

//~~Methods~~
string get_line_fasta(ifstream* partToExamine);
vector<string> parse_this_fasta(string fastaFile); //??? or void parse_this_fasta(string fastaFile,void (*insertionFunction)(string));
}
#endif
