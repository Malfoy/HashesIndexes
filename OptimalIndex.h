#ifndef OptimalIndex
#define optimalIndex

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <bit>
#include <bitset>
#include <stdint.h>
#include <utility>



using namespace std;



typedef  vector<uint32_t>  genome_list;



class Bucket{
public:
// Attributes
uint8_t maximal_fingerprint;
genome_list* fingerprint_to_list;

//~~Constructor~~
Bucket(uint8_t maximal_fingerprint){
    fingerprint_to_list=new genome_list[maximal_fingerprint];
}

~Bucket(){
    delete fingerprint_to_list;
}

//Methods
genome_list get_genomes(uint fp){
    return fingerprint_to_list[fp];
}

void add_genome(uint32_t genome_id,uint fp){
    fingerprint_to_list[fp].push_back(genome_id);
}
};



class OptimalIndex {
public:

//~~Attributes~~
uint64_t nb_minimizer;
uint size_fingerprint;
uint maximal_fingerprint;
vector<Bucket> minimizer_to_bucket; 
uint64_t nb_genomes;

//~~Constructor~~
OptimalIndex(int64 bucketsNumber,uint size_fingerprint = 8);

//~~Methods~~
void add_sequence(const string& theSequence);
vector<double> query_sequence(const string&str);
vector<pair<uint32_t,double>> query_sequence(const string&str,double threshold);

};

#endif