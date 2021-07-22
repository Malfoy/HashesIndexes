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



//~~Constructor~~
TestTable::TestTable(uint32 genomeQuantityForTest) : kmerSize(31), nbGenomes(genomeQuantityForTest)
{
        unordered_map<string, vector<uint32> > hashTable; // will contains genomes references for every kmer
}


//~~Methods~~
void TestTable::parse_fasta_for_refTable(const string& fileName)
{
        ifstream theRead(fileName);
        while(not theRead.eof()) //put sequences string in genome vector while it's not End of File
        {
                record_sequence(get_line_fasta(&theRead),nbGenomes);
        }
        return;
}


//to find a number (genome here) in a vector
bool TestTable::ask_genomes_vector(vector<uint32> genomesVector, uint32 wantedGenome)
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
}



vector<double> TestTable::query_belonging_genome(const string& sequenceStr, double thresholdJaccard)
{
        vector<double> allScores(nbGenomes,0);//vector containing numbers of hit with the number genome for the index
        int sequenceSize(sequenceStr.size()), position(0), kmerSum(sequenceSize - kmerSize);
        for (position = 0; position < (kmerSum); position++) //comparison loop of kmer and select the best result.
        {
                string kmerToStudy{sequenceStr.substr (position,kmerSize)};
                if (hashTable.count(kmerToStudy) == 1)
                {
                        for (positionGen = 0; positionGen < hashTable[kmerToStudy].second.size(); positionGen++) //browsing of every genome associated with kmer
                        {
                                allScores[hashTable[kmerToStudy].second[positionGen]]++; //increment the value at the "genome number" index of vector allScores.
                        }
                }

        }
        assert(allScores[0] > 0); //to verify if there is the genome 0.

        for (positionGen=0; positionGen < nbGenomes; positionGen++) //browse the vector to transform hits number in Jaccard Index
        {
                allScores[positionGen] = allScores[positionGen]/kmerSum;
        }
        return allScores;
}
