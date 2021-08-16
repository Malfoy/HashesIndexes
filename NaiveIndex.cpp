#include <stdio.h>
#include <fstream>
#include <bitset>
#include<bits/stdc++.h>
#include <iostream>
#include <vector>
#include <atomic>  // to isolate the calcul thread
#include <functional>
#include <mutex>  // to orchestrate thread
#include <string>
#include <stdint.h>
#include <math.h>
#include <cassert>
#include "NaiveIndex.h"



using namespace std;


hash<std::string> kmer_hasher;




//~~Constructor~~

NaiveIndex::NaiveIndex(uint64 bucketsNumber,uint16 decimal_for_lsb,uint16 nbgenomes) : nb_minimizer(bucketsNumber), nb_genomes(nbgenomes), decimal_lsb(decimal_for_lsb), number_of_LSB(log2(decimal_for_lsb)), bit_to_keep_minimizer(log2(bucketsNumber)), kmerSize(31), matrix(bucketsNumber)
{
        matrix.resize(bucketsNumber);
}




//~~Methods~~

//  ~~Public~~

void NaiveIndex::index_sequences_from_fasta(const string& fileName)
{
        ifstream theRead(fileName);
        while(not theRead.eof()) //put sequences string in genome vector while it's not End of File
        {
                nb_genomes++;
                add_sketch(compute_sketch(get_line_fasta_for_naive(&theRead)));
        }
}


vector<long double> NaiveIndex::query_sequence(string sequenceSearched, long double acceptanceTreshold)
{
        assert(acceptanceTreshold >= 0 && acceptanceTreshold <= 1);//verify the value of acceptance treshold
        vector<long double> allScores(nb_genomes,0); //initialize the vector which will contains result
        uint matrixSize(matrix[0].size()), noHitBuckets(0), hitsCounter(0); //different initializing
        vector<uint8> vectorisedQuery(compute_sketch(sequenceSearched)); //compute sketch the sequence that we searched

        for (uint genomeY = 0; genomeY < (matrixSize); genomeY++)//browse the genome lines
        {
                for (uint bucketX = 0; bucketX < nb_minimizer; bucketX++)//browse the query and genome sketch
                {
                         if(vectorisedQuery[bucketX] != (decimal_lsb - 1)) //because we don't want record unsignificant bucket, by defaut 255
                        {
                                if((matrix[bucketX][genomeY]) == vectorisedQuery[bucketX])
                                {
                                        hitsCounter++;
                                }
                                else
                                {
                                        noHitBuckets++; //will participe at Jaccard index calcul
                                }
                        }
                }
                long double occurentJaccardIndex = (long double) hitsCounter/(hitsCounter+noHitBuckets);
                if(occurentJaccardIndex >= acceptanceTreshold) //record occurent score just above the treshold jaccard index
                {
                        allScores[genomeY] = occurentJaccardIndex;
                }
                hitsCounter = 0;//reset for a new loop
                noHitBuckets = 0;
                occurentJaccardIndex = 0;
        }
        return allScores;
}


vector<pair<long double,uint16>> NaiveIndex::sort_scores(vector<long double> allScoresVector)
{
  vector<pair<long double,uint16>> sortedScoresVector;
  for (uint genomeCursor = 0; genomeCursor < allScoresVector.size(); genomeCursor++)
  {
        sortedScoresVector.push_back(make_pair(allScoresVector[genomeCursor],genomeCursor));
  }
  sort (sortedScoresVector.rbegin(), sortedScoresVector.rend()); //rbegin (and rend) for descending else it would be begin
  return sortedScoresVector;
}


void NaiveIndex::show_sorted_scores(vector<pair<long double,uint16>> sortedScoresVector, uint howManyScoresToShow) // 0 mean all scores
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
    cout << sortedScoresVector[positionScore].first << "               " << sortedScoresVector[positionScore].second << endl;
  }
}





//  ~~Private~~

string NaiveIndex::get_line_fasta_for_naive(ifstream* partToExamine)
{
        string line,justTheSequence;
        getline(*partToExamine,line);
        char caracter=partToExamine->peek();
        while(caracter!='>' and caracter!=EOF) //avoid line with '>' and End Of File
        {
                getline(*partToExamine,line);
                line.erase(std::remove(line.begin(), line.end(), '\n'), line.end());
                justTheSequence+=line;
                caracter=partToExamine->peek();
        }
        return justTheSequence;
}


vector<uint8> NaiveIndex::compute_sketch(string sequenceStr)
{
        vector<uint8> result(nb_minimizer,(pow(2,number_of_LSB)-1)); // the resulting vector
        uint sequenceSize(sequenceStr.size()), position(0);

        /* a For loop to hash every kmer of a sequence, distribute them in
        different bucket with their MSB on one byte (255 buckets) and record
        some LSB number of the minimal hash least*/

        for (position = 0; position <= (sequenceSize - kmerSize); position++)
        {
                string kmerToHash{sequenceStr.substr (position,kmerSize)};
                int64 hashOfKmer(kmer_hasher(get_complement_or_not(kmerToHash)));
                if(result[get_bucket(hashOfKmer)] > get_hash(hashOfKmer))
                {
                        result[get_bucket(hashOfKmer)] = get_hash(hashOfKmer);
                }
        }
        return result;
}


void NaiveIndex::add_sketch(const vector<uint8>& sketchToAdd)
{


        int sketchSize(sketchToAdd.size());
        for (int position = 0; position < sketchSize; position++)
        {
                matrix[position].push_back(sketchToAdd[position]);
        }
}




string NaiveIndex::get_complement_or_not(string sequenceToComplement) // transform in uppercase
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


uint64 NaiveIndex::get_bucket(uint64 primaryHash) //giving minimizer
{
        return primaryHash>>(64-bit_to_keep_minimizer); //bit shift to know bucket, first bits, most significant bit (MSB)
}


uint8 NaiveIndex::get_hash(uint64 primaryHash) //record the LSB number that we choose
{
        primaryHash%=(1<<(64-bit_to_keep_minimizer));
        return get_minhash(primaryHash);
}




uint8 NaiveIndex::get_minhash(uint64 primaryHash)
{
        return primaryHash%(decimal_lsb - 1); //to obtain the hash we want, last bits, least significant bit (LSB)
}





//  ~~alternative~~

uint8 NaiveIndex::get_hyperloglog(uint64 primaryHash)
{
        return log2(primaryHash); //to obtain the hash we want, last bits, least significant bit (LSB)
}


uint64 NaiveIndex::xs(uint64 y) //Hash xshort method
{
        y^=(y<<13);
        y^=(y>>17);
        return (y^=(y<<15));
}


uint8 NaiveIndex::get_double_hyperloglog(uint64 primaryHash)
{
        uint8 result=min((uint8)log2(primaryHash),(uint8)15);
        result<<=4;
        uint8 hll2(log2(xs(primaryHash)));
        return result+=min(hll2,(uint8)15);

}


uint8 NaiveIndex::get_popcount(uint64 primaryHash)
{
        return __builtin_popcount(primaryHash);
}


uint8 NaiveIndex::get_hyper_minhashX(uint64 primaryHash, uint8 hyperMinhashNumber)
{//accepted values == 62,53,44,35,26
        uint8 result(0), bitShift(hyperMinhashNumber%10);
        assert((int(pow(2,(8-bitShift)))) % 4 == 0 && bitShift < 7);
        if (bitShift == 2)
        {
                result=get_hyperloglog(primaryHash);
        }
        else
        {
                result=min(get_hyperloglog(primaryHash),(uint8)(pow(2,(8-bitShift))-1));
        }
        result<<=bitShift;
        uint8 bitPower(pow(2,bitShift));
        result+=primaryHash%bitPower;
        return result;
}
