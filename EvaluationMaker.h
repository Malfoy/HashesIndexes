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
        double success = 0;
        double falsePositive = 0;
        double falseNegative = 0;
        double tooMuchKmer = 0;
        double notEnoughKmer = 0;
};



// class about a structure to evaluate results
class ComparisonMatrix {
public:

//~~Attribute~~
uint32 matrixHeight;
vector<vector<double> > testMatrix;
vector<Comparison_scores> final_comparison_score_vector;

//~~Constructor~~
ComparisonMatrix();

//~~Methods~~
void add_result_vector(vector<double> resultVector);
void create_comparison();
void show_the_matrix();
};

#endif
