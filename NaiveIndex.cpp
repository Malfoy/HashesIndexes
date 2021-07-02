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



hash<std::string> kmer_hasher;


//constructor
NaiveIndex::NaiveIndex(int64 bucketsNumber,uint16 decimal_for_lsb) : nb_minimizer(bucketsNumber), nb_genomes(0), decimal_lsb(decimal_for_lsb), bit_to_keep_minimizer(log2(bucketsNumber)), kmerSize(31), matrix(bucketsNumber)
{
        matrix.resize(bucketsNumber);
}



//giving minimizer
uint64 NaiveIndex::get_bucket(int64 primaryHash)
{
        return primaryHash>>(64-bit_to_keep_minimizer); //bit shift to know bucket, first bits, most significant bit (MSB)
}



//record the LSB number that we choose
uint8 NaiveIndex::get_hash(int64 primaryHash)
{
        return primaryHash%decimal_lsb; //to obtain the hash we want, last bits, least significant bit (LSB)
}



vector<uint8> NaiveIndex::compute_sketch(const string& sequenceStr,const int kmerSize)
{
        vector<uint8> result(nb_minimizer,255); // the resulting vector

        // a For loop to hash every kmer of a sequence, distribute them in different
        // bucket with their MSB on one byte (255 buckets) and record some LSB number
        // of the minimal hash least

        int sequenceSize(sequenceStr.size()), position(0);

        for (position = 0; position < (sequenceSize - kmerSize); position++)
        {
                string kmerToHash{sequenceStr.substr (position,kmerSize)};
                int64 hashOfKmer(kmer_hasher(kmerToHash));
                if(result[get_bucket(hashOfKmer)] > get_hash(hashOfKmer))
                {
                        result[get_bucket(hashOfKmer)] =  get_hash(hashOfKmer);
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



vector<score_strct> NaiveIndex::query_sequence(const string& sequenceSearched, int acceptanceTreshold = 0)
{
        vector<score_strct> allScores;
        int matrixSize(matrix[0].size()), noHitBuckets(0), hitsCounter(0);
        vector<uint8> vectorisedQuery(compute_sketch(sequenceSearched, kmerSize));
        for (int genomeY = 0; genomeY < (matrixSize); genomeY++)//browse the buckets
        {
                for (int bucketX = 0; bucketX < nb_minimizer; bucketX++)//browse the query and genome sketch
                {
                        if(vectorisedQuery[bucketX] != 255) //because we don't want record unsignificant bucket
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
                if(hitsCounter > acceptanceTreshold)
                {
                        score_strct transitoryStructure;
                        transitoryStructure.genomeNumber = genomeY;
                        transitoryStructure.jaccardIndex = genomeY + (hitsCounter/((hitsCounter+noHitBuckets)*10)); // example : 101.029 => Bucket number 101 and 0.29 Jaccard Index;
                        allScores.push_back(transitoryStructure);
                }
                hitsCounter = 0;//reset for a new loop
        }
        return allScores;
}
