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
#include "NaiveIndex.h"


using namespace std;


//~~Constructor~~
ComparisonMatrix::ComparisonMatrix() : matrixHeight(0)
{
        vector<vector<double> > testMatrix;
        struct final_comparison_scores {
                vector<uint32> studiedGenomeNumber;
                vector<double> success;
                vector<double> falsePositive;
                vector<double> falseNegative;
                vector<double> tooMuchKmer;
                vector<double> notEnoughKmer;
        };
}


//~~Methods~~
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
        testMatrix.push_back(resultVector)
}


//the index of the final_comparison_scores vector indicate the studied program version
//the index 0 is the "reference results vector" from testtable
void ComparisonMatrix::create_comparison()
{
        final_comparison_scores.studiedGenomeNumber = matrixHeight;

        for (methodLine = 1; methodLine < testMatrix.size(); methodLine++)
        {
                vector<uint32> differentScores(6);
                for (genomeIndex = 0; genomeIndex < matrixHeight; genomeIndex)
                {
                        if (testMatrix[0][genomeIndex] == testMatrix[methodLine])// if it's the same Jaccard index or both 0 it's a success
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
                final_comparison_scores.success.push_back(differentScores[1]/matrixHeight);
                final_comparison_scores.falsePositive.push_back(differentScores[2]/matrixHeight);
                final_comparison_scores.falseNegative.push_back(differentScores[3]/matrixHeight);
                final_comparison_scores.tooMuchKmer.push_back(differentScores[4]/matrixHeight);
                final_comparison_scores.notEnoughKmer.push_back(differentScores[5]/matrixHeight);
        }
}


void ComparisonMatrix::print_vector(vector<double> vectorToPrint)
{
        for (uint32 position = 0; position < (vectorToPrint.size()); position++)
        {
                cout << vect[p] << ' ';
        }
        cout << endl;
}


void ComparisonMatrix::show_the_matrix()
{
        cout << "Number of:" << endl << endl <<  "   matches   ";
        print_vector(final_comparison_scores.success);
        cout << "Number of:" << endl << endl <<  "   False positives   ";
        print_vector(final_comparison_scores.falsePositive);
        cout << "Number of:" << endl << endl <<  "   False negatives   ";
        print_vector(final_comparison_scores.falseNegative);
        cout << "Number of:" << endl << endl <<  "   with more kmer   ";
        print_vector(final_comparison_scores.tooMuchKmer);
        cout << "Number of:" << endl << endl <<  "   with less kmer   ";
        print_vector(final_comparison_scores.notEnoughKmer);
}
