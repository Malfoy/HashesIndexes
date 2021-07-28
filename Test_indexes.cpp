
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
        string theFasta("10Bact.fa");
        TestTable refTable(10);
        NaiveIndex firstIndex(10,8);
        refTable.parse_fasta_for_refTable(theFasta);
        firstIndex.add_fasta_for_naive(theFasta);
        string otherFasta("10Bact.fa");
        string oneGenome(firstIndex.get_line_fasta_for_naive(&otherFasta));
        ComparisonMatrix firstMatrix;
        vector<double> TestResultVector(TestTable.query_belonging_genome(oneGenome));
        vector<double> NaiveResultVector(firstIndex.query_sequence(oneGenome));
        firstMatrix.add_result_vector(TestResultVector);
        firstMatrix.add_result_vector(NaiveResultVector);
        firstMatrix.create_comparison();
        firstMatrix.show_the_matrix();
        
        return 0;
}
