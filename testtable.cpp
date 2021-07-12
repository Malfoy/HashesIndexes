#include <iostream>
#include <vector>
#include <string>  // for substr
#include <unordered_map>  // for Hashtable
#include <algorithm>  // for sorting
#include <utility> // for pair
#include "NaiveIndex.h"
#include <cassert> //for assertion



using namespace std;




hash<std::string> kmer_hasher;

//constructor
TestTable::TestTable(uint32 genomeQuantityForTest) : kmerSize(31), nbGenomes(genomeQuantityForTest)
{
  unordered_map<string, vector<uint32> > hashTable; // will contains genomes references for every kmer
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

/* There are two ways here to have genome number, first the index, then the second
value of the pair*/
vector<pair <uint32, uint32>> TestTable::query_belonging_genome(const string& sequenceStr)
{
  vector<pair <uint32, uint32>> allScores(nbGenomes,0);//vector of pairs containing numbers of hit and corresponding genome number
  int sequenceSize(sequenceStr.size()), position(0);
  for (position = 0; position < (sequenceSize - kmerSize); position++) //comparison loop of kmer and select the best result.
  {
    string kmerToStudy{sequenceStr.substr (position,kmerSize)};
    if (hashTable.count(kmerToStudy) == 1)
    {
      for (positionGen = 0; positionGen < hashTable[kmerToStudy].second.size(); positionGen++) //browsing of every genome associated with kmer
      {
        allScores[hashTable[kmerToStudy].second[positionGen]].first++; //increment the value at the "genome number" index of vector allScores.
        allScores[hashTable[kmerToStudy].second[positionGen]].second = allScores[hashTable[kmerToStudy].second[positionGen]]; // we'll certainly need the genome value in second position.
      }
    }

  }
  assert((allScores[0].second > 0)==1); //to verify there is the genome 0.
  return allScores;
}
