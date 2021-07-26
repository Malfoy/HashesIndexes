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
#include <cassert>
#include "EvaluationMaker.h"



using namespace std;


// VECTOR STRUCT ?
struct  Comparison_scores{
        uint32 studiedGenomeNumber = 0;
        double success = 0;
        double falsePositive = 0;
        double falseNegative = 0;
        double tooMuchKme = 0;
        double notEnoughKmer = 0;
};


//~~Constructor~~
ComparisonMatrix::ComparisonMatrix() : matrixHeight(0), testMatrix(), final_comparison_score_vector()
{
}


//~~Methods~~
// the first result vector must come from testtable
void ComparisonMatrix::add_result_vector(vector<double> resultVector)
{
        //Verify if vectors are the same size
        if (matrixHeight == 0)
        {
                matrixHeight = resultVector.size();
        }
        else
        {
                assert(matrixHeight==resultVector.size());
        }
        testMatrix.push_back(resultVector);
}


//the index of the final_comparison_scores vector (methodLine) indicates the studied program version
//the index 0 is the "reference results vector" from testtable
void ComparisonMatrix::create_comparison(uint32 methodLine)
{
        final_comparison_score_vector[methodLine].studiedGenomeNumber = matrixHeight;
        long unsigned int methodLine(0);//this type because it's compared to size()
        uint32 genomeIndex(0);

        for (methodLine = 1; methodLine < testMatrix.size(); methodLine++)
        {
                vector<uint32> differentScores(6);//sum the corresponding result at his category (different score[1] = success score)
                for (genomeIndex = 0; genomeIndex < matrixHeight; genomeIndex++)
                {
                        if (testMatrix[0][genomeIndex] == testMatrix[methodLine][genomeIndex])// if it's the same Jaccard index or both 0 it's a success
                        {
                                differentScores[1]++;
                        }
                        else if (testMatrix[methodLine][genomeIndex] == 0)//false positive
                        {
                                differentScores[2]++;
                        }
                        else if (testMatrix[0][genomeIndex] == 0)//false negative
                        {
                                differentScores[3]++;
                        }
                        else if (testMatrix[methodLine][genomeIndex] > testMatrix[0][genomeIndex])//too much genome
                        {
                                differentScores[4]++;
                        }
                        else if (testMatrix[methodLine][genomeIndex] < testMatrix[0][genomeIndex])//not enough genome
                        {
                                differentScores[5]++;
                        }
                }
                final_comparison_score_vector[methodLine].success = (differentScores[1]/matrixHeight);
                final_comparison_score_vector[methodLine].falsePositive = (differentScores[2]/matrixHeight);
                final_comparison_score_vector[methodLine].falseNegative = (differentScores[3]/matrixHeight);
                final_comparison_score_vector[methodLine].tooMuchKmer = (differentScores[4]/matrixHeight);
                final_comparison_score_vector[methodLine].notEnoughKmer = (differentScores[5]/matrixHeight);
        }
}


void ComparisonMatrix::show_the_matrix()
{
  long unsigned int methodLine(0);//this type because it's compared to size()
  for (methodLine = 1; methodLine < testMatrix.size(); methodLine++)
  {
        cout << "For the method " << methodLine << endl;
        cout << "Number of:" << endl << endl <<  "   matches   ";
        cout << final_comparison_score_vector[methodLine].success << endl;
        cout << endl << endl <<  "   False positives   ";
        cout << final_comparison_score_vector[methodLine].falsePositive << endl;
        cout << endl << endl <<  "   False negatives   ";
        cout << final_comparison_score_vector[methodLine].falseNegative) << endl;
        cout << endl << endl <<  "   with more kmer   ";
        cout << final_comparison_score_vector[methodLine].tooMuchKmer) << endl;
        cout << endl << endl <<  "   with less kmer   ";
        cout << final_comparison_score_vector[methodLine].notEnoughKmer) << endl;
  }
}
