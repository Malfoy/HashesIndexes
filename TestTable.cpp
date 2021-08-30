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
TestTable::TestTable(uint32 genomeQuantityForTest) : goodGenomeVector({-2}), kmerSize(31), nbGenomes(genomeQuantityForTest), hashTable(), kmerCountVector()
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


vector<vector<long double>> TestTable::get_all_queries_scores_for_table(vector<vector<long double>> ultimVector, const string& fileName)
{
  ifstream theRead(fileName);
  vector<int> wichGoodGenome(goodGenomeVector);
  while(not theRead.eof()) //put sequences string in genome vector while it's not End of File
  {
          vector<pair<long double,uint16>> sortedVectorTest(sort_scores(query_belonging_genome(get_line_fasta_for_testtable(&theRead))));
          vector <long double> goodVector;
          if (wichGoodGenome[0] != -1)
          {
            if (wichGoodGenome[0] == -2)
            {
              wichGoodGenome[0] = sortedVectorTest[0].second;
              if (sortedVectorTest[0].first == 0)
              {
                wichGoodGenome[0] = 0;
              }
            }
            else
            {
              wichGoodGenome.push_back(sortedVectorTest[0].second);
              if (sortedVectorTest[0].first == 0)
              {
                wichGoodGenome.back() = 0;
              }
            }
          }
          for(uint position = 0; position < (uint) sortedVectorTest.size(); position++)
          {
            goodVector.push_back(sortedVectorTest[position].first);
          }
          ultimVector.push_back(goodVector);
  }
  goodGenomeVector = wichGoodGenome;
  return ultimVector;
}


vector<int> TestTable::get_good_genome_vector()
{
  return goodGenomeVector;
}


vector<long double> TestTable::query_belonging_genome(string sequenceStr)
{
        vector<long double> allScores(nbGenomes,0);//vector containing numbers of hit with the number genome for the index
        int sequenceSize(sequenceStr.size()), position(0), kmerhittoremove(0), kmerSumQuery(sequenceSize - kmerSize + 1);
        long unsigned int positionGen(0);//this type because it's compared to size()
        for (position = 0; position < (kmerSumQuery); position++) //comparison loop of kmer and select the best result.
        {
                string kmerToStudybeforecomplement{sequenceStr.substr (position,kmerSize)};
                string kmerToStudy(get_complement_or_not(kmerToStudybeforecomplement));
                if (hashTable.count(kmerToStudy) == 1)
                {
                        for (positionGen = 0; positionGen < hashTable[kmerToStudy].size(); positionGen++) //browsing of every genome associated with kmer
                        {
                                if (hashTable[kmerToStudy][positionGen].second == true)
                                {
                                  allScores[hashTable[kmerToStudy][positionGen].first]++; //increment the value at the "genome number" index of vector allScores.
                                  hashTable[kmerToStudy][positionGen].second = false; // prevent from counting gain the same kmer for the same genome
                                }
                                else
                                {
                                  kmerhittoremove++; // for the Jaccard index
                                }
                        }
                }
        }
        assert(allScores[0] > 0); //to verify if there is the genome 0.
        uint32 positionGeno(0);
        for (positionGeno=0; positionGeno < nbGenomes; positionGeno++) //browse the vector to transform hits number in Jaccard Index
        {
          allScores[positionGeno] = (long double) allScores[positionGeno]/(allScores[positionGeno] + (kmerSumQuery - kmerhittoremove - allScores[positionGeno]) + (kmerCountVector[positionGeno] - allScores[positionGeno]));
        }
        return allScores;
}


vector<pair<long double,uint16>> TestTable::sort_scores(vector<long double> allScoresVector, long double thresholdJaccard)
{
  assert(thresholdJaccard >= 0 && thresholdJaccard <= 1);//verify the value of acceptance treshold
  vector<pair<long double,uint16>> sortedScoresVector;
  uint thresholdPosition(0);
  for (uint genomeCursor = 0; genomeCursor < allScoresVector.size(); genomeCursor++)
  {
        sortedScoresVector.push_back(make_pair(allScoresVector[genomeCursor],genomeCursor));
  }
  sort (sortedScoresVector.rbegin(), sortedScoresVector.rend()); //rbegin (and rend) for descending else it would be begin
  for (thresholdPosition = 0; thresholdPosition < sortedScoresVector.size(); thresholdPosition++) // loop to keep juste genome with value abose treshold jaccard index
  {
    if (thresholdJaccard > sortedScoresVector[thresholdPosition].first) // if below treshold stop and record subvector with value above the treshold
    {
      vector<pair<long double,uint16>> trimSortedScoresVector{sortedScoresVector.begin(), sortedScoresVector.begin() + thresholdPosition};
      sortedScoresVector = trimSortedScoresVector;
      break;
    }
  }
  if (sortedScoresVector.size()==0)
  {
    cout << "   SortedScoresVector is empty, tresholdJaccard may be too high or no genome recorded. " << endl;
  }
  return sortedScoresVector;
}


void TestTable::show_sorted_scores(vector<pair<long double,uint16>> sortedScoresVector, uint howManyScoresToShow) // 0 mean all scores
{
  cout << endl << "  ~    TESTABLE SORTED SCORES WITH HIS GENOME NUMBER     ~  " << endl;
  cout << "Number Genome"<< "     " << "Jaccard Index" << endl;
  uint limitNumber(sortedScoresVector.size());
  if (limitNumber > howManyScoresToShow && howManyScoresToShow != 0)
  {
    limitNumber = howManyScoresToShow;
  }
  for (uint positionScore = 0; positionScore < limitNumber; positionScore++)
  {
    cout << sortedScoresVector[positionScore].second << "                    " << sortedScoresVector[positionScore].first << endl;
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
                line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                justTheSequence+=line;
                caracter=partToExamine->peek();
        }
        return justTheSequence;
}


void TestTable::record_sequence(string sequenceStr, const uint32 Genome)
{
        // insertion loop of kmer
        int sequenceSize(sequenceStr.size()), position(0), Kmercount(0);
        for (position = 0; position <= (sequenceSize - kmerSize); position++)
        {
                string kmerToStudybeforecomplement{sequenceStr.substr(position,kmerSize)};
                string kmerToStudy(get_complement_or_not(kmerToStudybeforecomplement));
                if (hashTable.count(kmerToStudy) == 1) // verify if there is already the kmer
                {
                        if (ask_genomes_vector(hashTable[kmerToStudy],Genome) == false) // if it's not the same last genome then it push back.
                        {
                                hashTable[kmerToStudy].push_back(make_pair(Genome, true));
                        }
                }
                else
                {
                        //hashTable[kmerToStudy].keptKmer = kmerToStudy;
                        hashTable[kmerToStudy].push_back(make_pair(Genome, true));
                        Kmercount++; // for the Jaccard index
                }
        }
        kmerCountVector.push_back(Kmercount); // for the Jaccard index
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




bool TestTable::ask_genomes_vector(vector<pair<uint32,bool>> genomesVector, uint32 wantedGenome)
{
        if (genomesVector.back().first == wantedGenome)
        {
                return true;
        }
        else
        {
                return false;
        }
}
