#ifndef TESTTABLE
#define TESTTABLE
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



class TestTable {
public:

//~~Attributes~~
uint32 Kmercount;
int64 kmerSize;
uint32 nbGenomes;
std::tr1::unordered_map<string, vector<pair<uint32,bool>>> hashTable;

//~~Constructor~~
TestTable(uint32 genomeQuantityForTest = 0);

//~~Methods~~
//    ~~public~~
void parse_fasta_for_refTable(const string& fileName);
vector<long double> query_belonging_genome(string sequenceStr, long double thresholdJaccard = 0);
vector<pair<long double,uint16>> sort_scores(vector<long double> allScoresVector);
void show_sorted_scores(vector<pair<long double,uint16>> sortedScoresVector, uint howManyScoresToShow = 0);

//    ~~private~~
string get_line_fasta_for_testtable(ifstream* partToExamine);
void record_sequence(string sequenceStr, const uint32 Genome);
string get_complement_or_not(string sequenceToComplement);

  bool ask_genomes_vector(vector<pair<uint32,bool>> genomesVector, uint32 wantedGenome); //to find a number (genome here) in a vector
};

#endif
