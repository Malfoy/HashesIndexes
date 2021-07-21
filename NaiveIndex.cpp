#include <stdio.h>
#include <fstream>
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
NaiveIndex::NaiveIndex(int64 bucketsNumber,uint16 decimal_for_lsb) : nb_minimizer(bucketsNumber), nb_genomes(0), decimal_lsb(decimal_for_lsb), bit_to_keep_minimizer(log2(bucketsNumber)), kmerSize(31), matrix(bucketsNumber)
{
        matrix.resize(bucketsNumber);
}


//~~Methods~~
//giving minimizer
uint64 NaiveIndex::get_bucket(int64 primaryHash)
{
        return primaryHash>>(64-bit_to_keep_minimizer); //bit shift to know bucket, first bits, most significant bit (MSB)
}



uint8 NaiveIndex::get_hash(int64 primaryHash)
{
        primaryHash%=(1<<(64-bit_to_keep_minimizer));
        return get_minhash(primaryHash);
}



//record the LSB number that we choose
uint8 NaiveIndex::get_minhash(int64 primaryHash)
{
        return primaryHash%decimal_lsb; //to obtain the hash we want, last bits, least significant bit (LSB)
}


uint8 NaiveIndex::get_hyperloglog(int64 primaryHash)
{
        return log2(primaryHash); //to obtain the hash we want, last bits, least significant bit (LSB)
}


//Hash xshort method
uint64_t xs(uint64_t y){
        y^=(y<<13);
        y^=(y>>17);
        return (y^=(y<<15));
}



uint8 NaiveIndex::get_double_hyperloglog(int64 primaryHash)
{
        uint8 result=min((uint8)log2(primaryHash),(uint8)15);
        result<<=4;
        uint8 hll2(log2(xs(primaryHash)));
        return result+=min(hll2,(uint8)15);

}


uint8 NaiveIndex::get_popcount(int64 primaryHash)
{
        return __builtin_popcount(primaryHash);
}


uint8 NaiveIndex::get_hyper_minhashX(int64 primaryHash, uint8 hyperMinhashNumber)
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


vector<double> NaiveIndex::query_sequence(const string& sequenceSearched, double acceptanceTreshold)
{
        assert(acceptanceTreshold >= 0 && acceptanceTreshold <= 1);
        vector<double> allScores(nb_genomes,0);
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
                double occurentJaccardIndex = hitsCounter/(hitsCounter+noHitBuckets);
                if(occurentJaccardIndex > acceptanceTreshold) //record occurent score just above the treshold jaccard index
                {
                        allScores[genomeY] = occurentJaccardIndex;
                }
                hitsCounter = 0;//reset for a new loop
                occurentJaccardIndex = 0;
        }
        return allScores;
}
