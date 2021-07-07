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


void TestTable::record_sequence(const string& sequenceStr, const string& Genome)
{
  //insertion loop of kmer
  int sequenceSize(sequenceStr.size()), position(0);

  for (position = 0; position < (sequenceSize - kmerSize); position++)
  {
          string kmerToStudy{sequenceStr.substr (position,kmerSize)};
          if hashTable[kmerToStudy].keptKmer == kmerToStudy && hashTable[kmerToStudy].genomeNumber //TODO mettre une liste
          {

          }

}

uint64 TestTable::query_belonging_genome(const string& sequenceStr)
{
  //comparison loop of kmer and select the best result.
  return genome_number;
}
