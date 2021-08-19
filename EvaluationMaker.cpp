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




//~~Constructor~~
ComparisonMatrix::ComparisonMatrix() : matrixHeight(0), testMatrix(), final_comparison_score_vector()
{
}


//~~Methods~~


//  ~~Public~~

// the first result vector must come from testtable
void ComparisonMatrix::add_result_vector(vector<long double> resultVector)
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
        final_comparison_score_vector.resize(testMatrix.size());
}


//the index of the final_comparison_scores vector (methodLine) indicates the studied program version
//the index 0 is the "reference results vector" from testtable
void ComparisonMatrix::create_comparison()
{
        long unsigned int methodLine(0);//this type because it's compared to size()
        final_comparison_score_vector[methodLine].studiedGenomeNumber = matrixHeight;
        uint32 genomeIndex(0);
        for (methodLine = 1; methodLine < testMatrix.size(); methodLine++)
        {
                vector<uint32> differentScores(7);//sum the corresponding result at his category (different score[1] = success score)
                for (genomeIndex = 0; genomeIndex < matrixHeight; genomeIndex++)
                {
                        if (testMatrix[0][genomeIndex] == testMatrix[methodLine][genomeIndex])// if it's the same Jaccard index or both 0 it's a success
                        {
                                differentScores[1]++;
                                if (testMatrix[0][genomeIndex] == 0)
                                {
                                  differentScores[6]++;
                                }
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
                final_comparison_score_vector[methodLine].success = (long double) differentScores[1]/matrixHeight;
                final_comparison_score_vector[methodLine].falsePositive = (long double) differentScores[2]/matrixHeight;
                final_comparison_score_vector[methodLine].falseNegative = (long double) differentScores[3]/matrixHeight;
                final_comparison_score_vector[methodLine].tooMuchKmer = (long double) differentScores[4]/matrixHeight;
                final_comparison_score_vector[methodLine].notEnoughKmer = (long double) differentScores[5]/matrixHeight;
                final_comparison_score_vector[methodLine].zeroHit = differentScores[6];

        }
}


void ComparisonMatrix::show_the_matrix()
{
  long unsigned int methodLine(0);//this type because it's compared to size()
  char space = ' ';
  for (methodLine = 1; methodLine < testMatrix.size(); methodLine++)
  {
        cout << endl << endl << endl << "Scores comparison with the Test table For the method " << methodLine << endl << "*************************************************" << endl;
        cout << "For " << matrixHeight << " genomes , percent of:" << endl << endl <<  "   Matches           ";
        cout << final_comparison_score_vector[methodLine].success * 100 << " %" << string(20 - get_number_digits(final_comparison_score_vector[methodLine].success * 100), space) << "    (" << final_comparison_score_vector[methodLine].success * matrixHeight << " genomes qualified)     (Full no hits: " << final_comparison_score_vector[methodLine].zeroHit << ")" << endl;
        cout << endl << endl <<  "   False positives   ";
        cout << final_comparison_score_vector[methodLine].falsePositive * 100 << " %" << string(20 - get_number_digits(final_comparison_score_vector[methodLine].falsePositive * 100), space) << "    (" << final_comparison_score_vector[methodLine].falsePositive * matrixHeight << " genomes qualified)" <<  endl;
        cout << endl << endl <<  "   False negatives   ";
        cout << final_comparison_score_vector[methodLine].falseNegative * 100 << " %" << string(20 - get_number_digits(final_comparison_score_vector[methodLine].falseNegative * 100), space) << "    (" << final_comparison_score_vector[methodLine].falseNegative * matrixHeight << " genomes qualified)" <<  endl;
        cout << endl << endl <<  "   With more kmer    ";
        cout << final_comparison_score_vector[methodLine].tooMuchKmer * 100 << " %" << string(20 - get_number_digits(final_comparison_score_vector[methodLine].tooMuchKmer * 100), space) << "    (" << final_comparison_score_vector[methodLine].tooMuchKmer * matrixHeight << " genomes qualified)" <<  endl;
        cout << endl << endl <<  "   With less kmer    ";
        cout << final_comparison_score_vector[methodLine].notEnoughKmer * 100 << " %" << string(20 - get_number_digits(final_comparison_score_vector[methodLine].notEnoughKmer * 100), space) << "    (" << final_comparison_score_vector[methodLine].notEnoughKmer * matrixHeight << " genomes qualified)" <<  endl << endl << "*************************************************" << endl << endl;
  }
}


void ComparisonMatrix::write_result(string fileName)
{
  ofstream myfile;
  myfile.open (fileName);
  myfile << "Method_Number,Success,False_positive,False_negative,with_too_much_kmer,with_not_enough_kmer" << endl;
  for (uint methodPosition = 0; methodPosition < (uint) final_comparison_score_vector.size(); methodPosition++)
  {
    myfile << methodPosition << "," << final_comparison_score_vector[methodPosition].success << "," << final_comparison_score_vector[methodPosition].falsePositive << "," << final_comparison_score_vector[methodPosition].falseNegative << "," << final_comparison_score_vector[methodPosition].tooMuchKmer << "," << final_comparison_score_vector[methodPosition].notEnoughKmer << endl;
  }
  myfile.close();
}





//  ~~Private~~


int ComparisonMatrix::get_number_digits(long double aNumber)//use for number spaces for show the matrix
{
    int digits(0), limit(0);
    while ( (int) aNumber != 0)
    {
      aNumber /= 10;
      digits++;
    }
    while (aNumber - (int) aNumber && limit != 5 && aNumber != 0)
    {
      aNumber *=10;
      digits++;
      limit++;
    }
    if (digits==0)
    {
      return 1;
    }
    else
    {
      return digits;
    }
}
