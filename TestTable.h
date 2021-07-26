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



class TestTable {
public:

//~~Attributes~~
int64 kmerSize;
uint32 nbGenomes;
std::tr1::unordered_map<string, vector<uint32> > hashTable;

//~~Constructor~~
TestTable(uint32 genomeQuantityForTest);

//~~Methods~~
bool ask_genomes_vector(vector<uint32> genomesVector, uint32 wantedGenome);
void record_sequence(const string& sequenceStr, const uint32 Genome);
void parse_fasta_for_refTable(const string& fileName);
vector<double> query_belonging_genome(const string& sequenceStr, double thresholdJaccard);
string get_line_fasta(ifstream* partToExamine);
};

#endif