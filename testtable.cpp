#include <iostream>
#include <string>
#include <unordered_map>
#include "NaiveIndex.h"

using namespace std;

hash<std::string> kmer_hasher;

//constructor
TestTable::TestTable() : kmerSize(31)
{
  unordered_map<string, hashTable_cell_strct > hashTable;
}


//to find a number (genome here) in a vector
//method found on https://thispointer.com/c-how-to-find-an-element-in-vector-and-get-its-index/
bool testtable::ask_genomes_vector(vector<uint32> genomesVector, uint32 wantedGenome)
{
  vector<uint32>::iterator genIter = find_if(genomesVector.begin(), genomesVector.end(), [](const uint32 & val){
                                                                                            if (val == wantedGenome)
                                                                                                return true;
                                                                                            return false;
                                                                                        });
}


void TestTable::record_sequence(const string& sequenceStr, const uint32 Genome)
{
  //insertion loop of kmer
  int sequenceSize(sequenceStr.size()), position(0);

  for (position = 0; position < (sequenceSize - kmerSize); position++)
  {
          string kmerToStudy{sequenceStr.substr (position,kmerSize)};
          if (hashTable[kmerToStudy].keptKmer == kmerToStudy)
          {
            if (ask_genomes_vector(hashTable[kmerToStudy].genomesNumberV,Genome) == false)
            {
              hashTable[kmerToStudy].genomesNumberV.push_back(Genome);
            }
          }
          else
          {
            hashTable[kmerToStudy].keptKmer = kmerToStudy;
            hashTable[kmerToStudy].genomesNumberV.push_back(Genome);
          }
}

uint64 TestTable::query_belonging_genome(const string& sequenceStr)
{
  //comparison loop of kmer and select the best result.
  return genome_number;
}
