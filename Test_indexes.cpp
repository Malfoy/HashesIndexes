
#include <algorithm>
#include <array>
#include <chrono>
#include <cstring>
#include <string>
#include <fstream>
#include <functional>
#include <iostream>
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "NaiveIndex.h"
#include "EvaluationMaker.h"
#include "TestTable.h"






int main(int argc, char** argv)
{
        string theFasta("fastaest.fa");
        cout << endl << "the fasta" << endl;
        TestTable refTable(2);//Don't forget in naive index genome number
        cout << endl << "reftable" << endl;
        NaiveIndex firstIndex(2,8);
        cout << endl << "firstindex" << endl;
        refTable.parse_fasta_for_refTable(theFasta);
        cout << endl << "parse_fasta_for_refTable" << endl;
        firstIndex.add_fasta_for_naive(theFasta);
        cout << endl << "add_fasta_for_naive" << endl;
        string otherFasta("fastest.fa");
        cout << endl << "otherFasta" << endl;
        ifstream oneRead(otherFasta);
        cout << endl << "oneRead" << endl;
        string oneGenome("azertyuiopqsdfghjklmwxcvbnbvcxwmlkjhgfdsqpoiuytreza");
        cout << endl << "oneGenome extracted" << endl;
        ComparisonMatrix firstMatrix;
        cout << endl << "firstMatrix" << endl;
        vector<double> TestResultVector(refTable.query_belonging_genome(oneGenome));
        cout << endl << "TestResultVector" << endl;
        vector<double> NaiveResultVector(firstIndex.query_sequence(oneGenome));
        cout << endl << "NaiveResultVector" << endl;
        firstMatrix.add_result_vector(TestResultVector);
        cout << endl << "firstMatrix.add_result_vector(TestResultVector)" << endl;
        firstMatrix.add_result_vector(NaiveResultVector);
        cout << endl << "firstMatrix.add_result_vector(NaiveResultVector)" << endl;
        firstMatrix.create_comparison();
        cout << endl << "firstMatrix.create_comparison()" << endl;
        firstMatrix.show_the_matrix();
        cout << endl << "firstMatrix.show_the_matrix()" << endl;

        return 0;
}
