#ifndef COMPARISON
#define COMPARISON
#define NDEBUG  //ASSERT ON TESTTABLE.CPP

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
#include <tr1/unordered_map>
#include <chrono>

// removal '_t'
typedef int8_t int8;
typedef uint8_t uint8;
typedef int16_t int16;
typedef uint16_t uint16;
typedef int32_t int32;
typedef uint32_t uint32;
typedef int64_t int64;
typedef uint64_t uint64;



using namespace std;



struct  Comparison_scores {
        uint32 studiedGenomeNumber = 0;
        long double success = 0;
        long double falsePositive = 0;
        long double falseNegative = 0;
        long double tooMuchKmer = 0;
        long double notEnoughKmer = 0;
        uint32 zeroHit = 0;
        double timeForDataBase = 999;
        double timeForQuery = 999;
};



// class about a structure to evaluate results
class ComparisonMatrix {
public:

//~~Attribute~~
uint32 matrixHeight;
vector<vector<long double> > testMatrix;
vector<Comparison_scores> final_comparison_score_vector;
//~~Constructor~~
ComparisonMatrix();

//~~Methods~~
//    ~~public~~
void add_result_vector(vector<long double> resultVector);
void create_comparison();
void show_the_matrix();
void write_result(string fileName, string recordChoice = "comparison", bool timeOption = 0); // "comparison", "jaccard", "jaccardBucketsNumberdifference", "jaccardBucketsNumberdifference(col_inversion)";
void fill_time(uint32 comparisonToFill, std::chrono::duration<double> TimeDatabase, std::chrono::duration<double> TimeQuery);
string get_string_two_best_JacInd(); //use for Error rate list experiment
string get_other_string_from_ultim_vector(vector<vector<long double>> ultimVector, uint queriesNumber,uint cursor); //use for thousands genome experience
void compare_good_genome(vector<int> refVector, vector<int> indVector); //use for thousands genome experience

//    ~~private~~
int get_number_digits(long double aNumber);
};

#endif
