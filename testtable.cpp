#include <iostream>
#include <string>
#include <unordered_map>
#include "NaiveIndex.h"

using namespace std;

hash<std::string> kmer_hasher;

//constructor
TestTable::TestTable() : kmerSize(31)
{
  unordered_map<std::string, std::int> hashTable;
}


void TestTable::record_sequence(const string& sequenceStr)
{
  //insertion loop of kmer
}

uint64 TestTable::query_belonging_genome(const string& sequenceStr)
{
  //comparison loop of kmer and select the best result.
  return genome_number;
}
