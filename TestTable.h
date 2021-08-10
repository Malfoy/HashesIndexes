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
int64 kmerSize;
uint32 nbGenomes;
std::tr1::unordered_map<string, vector<uint32> > hashTable;

//~~Constructor~~
TestTable(uint32 genomeQuantityForTest);

//~~Methods~~
//    ~~public~~
void parse_fasta_for_refTable(const string& fileName);
vector<double> query_belonging_genome(const string& sequenceStr, double thresholdJaccard = 0);

//    ~~private~~
string get_line_fasta_for_testtable(ifstream* partToExamine);
void record_sequence(const string& sequenceStr, const uint32 Genome);

  bool ask_genomes_vector(vector<uint32> genomesVector, uint32 wantedGenome); //to find a number (genome here) in a vector
};

#endif
