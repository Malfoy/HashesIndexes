#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <functional>
#include <mutex>
#include <string>
#include <stdint.h>
#include <math.h>
#include "NaiveIndex.h"



using namespace std;



hash<std::string> kmerHasher;



NaiveIndex::NaiveIndex(int64_t bucketsNumber) : nb_minimizer(bucketsNumber), nb_genomes(0), bit_to_keep_minimizer(log2(bucketsNumber)), kmerSize(31), matrix(bucketsNumber)
{
        matrix.resize(bucketsNumber);
}



uint64_t NaiveIndex::get_bucket(int64_t primaryHash)
{
        return primaryHash>>(64-bit_to_keep_minimizer); //bit shift to know bucket, first bits, most significant bit (MSB)
}



uint8_t NaiveIndex::get_hash(int64_t primaryHash)
{
        return primaryHash%256; //to obtain the hash we want, last bits, least significant bit (LSB)
}



vector<uint8_t> NaiveIndex::compute_sketch(const string& sequenceStr,const int kmerSize)
{
        vector<uint8_t> result(nb_minimizer,255); // the resulting vector

        // a For loop to hash every kmer of a sequence, distribute them in different
        // bucket with their MSB on one byte (255 buckets) and record some LSB number
        // of the minimal hash least

        int sequenceSize(sequenceStr.size()), position(0);

        for (position = 0; position < (sequenceSize - kmerSize); position++)
        {
                string kmerToHash{sequenceStr.substr (position,kmerSize)};
                int64_t hashOfKmer(kmerHasher(kmerToHash));
                if(result[get_bucket(hashOfKmer)] > get_hash(hashOfKmer))
                {
                        result[get_bucket(hashOfKmer)] =  get_hash(hashOfKmer);
                }
        }

        return result;
}

void NaiveIndex::add_sketch(const vector<uint8_t>& sketchToAdd)
{
        int sketchSize(sketchToAdd.size());
        for (int position = 0; position < sketchSize; position++)
        {
                matrix[position].push_back(sketchToAdd[position]);
        }
}

/* TO_ERASE_? bool compare_sketch(const vector<uint8_t>& sequenceStudied,vector<vector<uint8_t>> genomeNumber)
   {
   for (int position = 0 ; position < (sketchToAdd.size()) ; position++)
   {
    matrix[position].push_back(sketchToAdd[position]);
   }
   } */


vector<double> NaiveIndex::query_sequence(const string& sequenceSearched)
{
        int matrixSize(matrix[0].size());
        double hitsCounter(0);
        vector<double> bestScore(0);
        vector<uint8_t> vectorisedQuery(compute_sketch(sequenceSearched, kmerSize));//TODO
        for (int positionY = 0; positionY < (matrixSize); positionY++)//browse the buckets
        {
                for (int positionX = 0; positionX < nb_minimizer; positionX++)//browse the query
                {
                        if(vectorisedQuery[positionX] != 255 && (matrix[positionX][positionY]) == vectorisedQuery[positionX])
                        {
                                hitsCounter++;
                        }
                }
                if(hitsCounter > 0)//TODO maybe other trigger count
                {
                        bestScore.push_back(hitsCounter);
                }
                hitsCounter = 0;//reset the counter
        }
        return bestScore;
}
