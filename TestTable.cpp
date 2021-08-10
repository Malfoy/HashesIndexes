#include <iostream>
#include <vector>
#include <string>  // for substr
#include <unordered_map>  // for TestTable
#include <algorithm>  // for sorting
#include <utility> // for pair
#include "TestTable.h"
#include <cassert> //for assertion



using namespace std;








//~~Constructor~~
TestTable::TestTable(uint32 genomeQuantityForTest) : kmerSize(31), nbGenomes(genomeQuantityForTest), hashTable()
{
}


//~~Methods~~

//  ~~Public~~

void TestTable::parse_fasta_for_refTable(const string& fileName)
{
        ifstream theRead(fileName);
        while(not theRead.eof()) //put sequences string in genome vector while it's not End of File
        {
                record_sequence(get_line_fasta_for_testtable(&theRead),nbGenomes);
        }
        return;
}


vector<double> TestTable::query_belonging_genome(const string& sequenceStr, double thresholdJaccard)
{
        vector<double> allScores(nbGenomes,0);//vector containing numbers of hit with the number genome for the index
        int sequenceSize(sequenceStr.size()), position(0), kmerSum(sequenceSize - kmerSize);
        long unsigned int positionGen(0);//this type because it's compared to size()
        for (position = 0; position < (kmerSum); position++) //comparison loop of kmer and select the best result.
        {
                string kmerToStudy{sequenceStr.substr (position,kmerSize)};
                if (hashTable.count(kmerToStudy) == 1)
                {
                        for (positionGen = 0; positionGen < hashTable[kmerToStudy].size(); positionGen++) //browsing of every genome associated with kmer
                        {
                                allScores[hashTable[kmerToStudy][positionGen]]++; //increment the value at the "genome number" index of vector allScores.
                        }
                }

        }
        assert(allScores[0] > 0); //to verify if there is the genome 0.
        uint32 positionGeno(0);
        for (positionGeno=0; positionGeno < nbGenomes; positionGeno++) //browse the vector to transform hits number in Jaccard Index
        {
                allScores[positionGen] = allScores[positionGen]/kmerSum;
        }
        return allScores;
}





//  ~~Private~~

string TestTable::get_line_fasta_for_testtable(ifstream* partToExamine)
{
        string line,justTheSequence;
        getline(*partToExamine,line);
        char caracter=partToExamine->peek();
        while(caracter!='>' and caracter!=EOF) //avoid line with '>' and End Of File
        {
                getline(*partToExamine,line);
                justTheSequence+=line;
                caracter=partToExamine->peek();
        }
        return justTheSequence;
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
                        if (ask_genomes_vector(hashTable[kmerToStudy],Genome) == false)
                        {
                                hashTable[kmerToStudy].push_back(Genome);
                        }
                }
                else
                {
                        //hashTable[kmerToStudy].keptKmer = kmerToStudy;
                        hashTable[kmerToStudy].push_back(Genome);
                }
        }
}




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
