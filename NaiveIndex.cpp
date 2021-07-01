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



vector<double> NaiveIndex::query_sequence(const string& sequenceSearched, int acceptanceTreshold = 0)
{
        int matrixSize(matrix[0].size()), noHitBuckets(0), hitsCounter(0);
        double genomeNumberDotJaccarind(0);
        vector<double> bestScores(0);
        vector<uint8_t> vectorisedQuery(compute_sketch(sequenceSearched, kmerSize));
        for (int genomeY = 0; genomeY < (matrixSize); genomeY++)//browse the buckets
        {
                for (int bucketX = 0; bucketX < nb_minimizer; bucketX++)//browse the query
                {
                        if(vectorisedQuery[bucketX] != 255)
                        {
                                if((matrix[bucketX][genomeY]) == vectorisedQuery[bucketX])
                                {
                                        hitsCounter++;
                                }
                                else
                                {
                                        noHitBuckets++;
                                }
                        }
                }
                if(hitsCounter > acceptanceTreshold)
                {
                        genomeNumberDotJaccarind = genomeY + (hitsCounter/((hitsCounter+noHitBuckets)*10)); // example : 101.029 => Bucket number 101 and 0.29 Jaccard Index
                        bestScores.push_back(genomeNumberDotJaccarind);
                        genomeNumberDotJaccarind = 0;//reset for a new loop
                }
                hitsCounter = 0;//reset for a new loop
        }
        return bestScores;
}
