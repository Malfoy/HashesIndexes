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
                nbGenomes++;
        }
}


vector<long double> TestTable::query_belonging_genome(string sequenceStrbeforeComplement, long double thresholdJaccard)
{
        string sequenceStr(get_complement_or_not(sequenceStrbeforeComplement));
        vector<long double> allScores(nbGenomes,0);//vector containing numbers of hit with the number genome for the index
        int sequenceSize(sequenceStr.size()), position(0), kmerSum(sequenceSize - kmerSize + 1);
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
                allScores[positionGeno] = (long double) allScores[positionGeno]/(kmerSum);
        }
        return allScores;
}


vector<pair<long double,uint16>> TestTable::sort_scores(vector<long double> allScoresVector)
{
  vector<pair<long double,uint16>> sortedScoresVector;
  for (uint genomeCursor = 0; genomeCursor < allScoresVector.size(); genomeCursor++)
  {
        sortedScoresVector.push_back(make_pair(allScoresVector[genomeCursor],genomeCursor));
  }
  sort (sortedScoresVector.rbegin(), sortedScoresVector.rend()); //rbegin (and rend) for descending else it would be begin
  return sortedScoresVector;
}


void TestTable::show_sorted_scores(vector<pair<long double,uint16>> sortedScoresVector, uint howManyScoresToShow) // 0 mean all scores
{
  cout << "  ~    SORTED SCORES WITH HIS GENOME NUMBER     ~  " << endl;
  cout << "Jccrd Idx"<< "     " << "Genome number" << endl;
  uint limitNumber(sortedScoresVector.size());
  if (limitNumber > howManyScoresToShow && howManyScoresToShow != 0)
  {
    limitNumber = howManyScoresToShow;
  }
  for (uint positionScore = 0; positionScore < limitNumber; positionScore++)
  {
    cout << sortedScoresVector[positionScore].first << "               " << sortedScoresVector[positionScore].second << endl; //Why '-nan' without uint ??
  }
}





//  ~~Private~~

string TestTable::get_line_fasta_for_testtable(ifstream* partToExamine)
{
        string line,justTheSequence;
        getline(*partToExamine,line);
        char caracter=partToExamine->peek();
        while(caracter!='>' and caracter!=EOF) // avoid line with '>' and End Of File
        {
                getline(*partToExamine,line);
                justTheSequence+=line;
                caracter=partToExamine->peek();
        }
        return justTheSequence;
}


void TestTable::record_sequence(string sequenceStrbeforeComplement, const uint32 Genome)
{
        string sequenceStr(get_complement_or_not(sequenceStrbeforeComplement));
        // insertion loop of kmer
        int sequenceSize(sequenceStr.size()), position(0);
        for (position = 0; position < (sequenceSize - kmerSize); position++)
        {
                string kmerToStudy{sequenceStr.substr(position,kmerSize)};
                if (hashTable.count(kmerToStudy) == 1) // verify if there is already the kmer
                {
                        if (ask_genomes_vector(hashTable[kmerToStudy],Genome) == false) // if it's not the same last genome then it push back.
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


string TestTable::get_complement_or_not(string sequenceToComplement) // transform in uppercase
{
  string reverscomp("");
  for (int i=sequenceToComplement.length()-1; i>=0; i--)
  {
     switch (sequenceToComplement[i])
    {

      case 'A':
        reverscomp.append("T");
        break;
      case 'T':
        reverscomp.append("A");
        break;
      case 'G':
        reverscomp.append("C");
        break;
      case 'C':
        reverscomp.append("G");
        break;
      case 'a':
        replace(sequenceToComplement.begin(),sequenceToComplement.end(),'a','A');
        reverscomp.append("T");
        break;
      case 'c':
        replace(sequenceToComplement.begin(),sequenceToComplement.end(),'c','C');
        reverscomp.append("G");
        break;
      case 'g':
        replace(sequenceToComplement.begin(),sequenceToComplement.end(),'g','G');
        reverscomp.append("C");
        break;
      case 't':
        replace(sequenceToComplement.begin(),sequenceToComplement.end(),'t','T');
        reverscomp.append("A");
        break;
    }
  }

  if(sequenceToComplement < reverscomp)
  {
    return sequenceToComplement;
  }
  else
  {
    return reverscomp;
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
