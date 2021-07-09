#include <iostream>
#include <string>
#include <unordered_map>
#include "NaiveIndex.h"

using namespace std;

hash<std::string> kmer_hasher;

//constructor
TestTable::TestTable() : kmerSize(31)
{
  unordered_map<string, vector<uint32> > hashTable;
}


//to find a number (genome here) in a vector
bool Testtable::ask_genomes_vector(vector<uint32> genomesVector, uint32 wantedGenome)
{
  if (genomesVector.back() == wantedGenome)
  {
    return true;
  }
  else
  {
    return false;
  }
}


void TestTable::record_sequence(const string& sequenceStr, const uint32 Genome)
{
  //insertion loop of kmer
  int sequenceSize(sequenceStr.size()), position(0);
  for (position = 0; position < (sequenceSize - kmerSize); position++)
  {
          string kmerToStudy{sequenceStr.substr (position,kmerSize)};
          if (hashTable.count(kmerToStudy) == 1)
          {
            if (ask_genomes_vector(hashTable[kmerToStudy].second,Genome) == false)
            {
              hashTable[kmerToStudy].second.push_back(Genome);
            }
          }
          else
          {
            //hashTable[kmerToStudy].keptKmer = kmerToStudy;
            hashTable[kmerToStudy].second.push_back(Genome);
          }
}

uint64 TestTable::query_belonging_genome(const string& sequenceStr)
{
  //comparison loop of kmer and select the best result.
  return genome_number;
}
