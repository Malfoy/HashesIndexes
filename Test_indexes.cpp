
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
        //time
        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        string theFasta("plouf.fa");
        cout << endl << "the fasta" << endl;
        TestTable refTable;
        cout << endl << "reftable" << endl;
        NaiveIndex firstIndex(524288,256);//((bucketnumbers, binary decimal (2,4,8...1024...) for record))
        cout << endl << "firstindex" << endl;
        refTable.parse_fasta_for_refTable(theFasta);
        cout << endl << "parse_fasta_for_refTable" << endl;
        firstIndex.index_sequences_from_fasta(theFasta);
        cout << endl << "add_fasta_for_naive" << endl;
        ifstream theFileQuery("tada.fa");
        string oneGenome(firstIndex.get_line_fasta_for_naive(&theFileQuery));
        cout << endl << "tada extrait" << endl;
        ComparisonMatrix firstMatrix;
        cout << endl << "firstMatrix" << endl;
        vector<long double> TestResultVector(refTable.query_belonging_genome(oneGenome));
        cout << endl << "TestResultVector" << endl;
        vector<pair<long double,uint16>> sortedVectorTest(refTable.sort_scores(TestResultVector,0.1));
        cout << endl << "sort vector test" << endl;
        refTable.show_sorted_scores(sortedVectorTest);
        vector<long double> NaiveResultVector(firstIndex.query_sequence(oneGenome));
        cout << endl << "NaiveResultVector" << endl;
        vector<pair<long double,uint16>> sortedVector(firstIndex.sort_scores(NaiveResultVector,0.1));
        cout << endl << "sort vector" << endl;
        firstIndex.show_sorted_scores(sortedVector);
        cout << endl << "show sorted vector" << endl;
        firstMatrix.add_result_vector(TestResultVector);
        cout << endl << "firstMatrix.add_result_vector(TestResultVector)" << endl;
        firstMatrix.add_result_vector(NaiveResultVector);
        cout << endl << "firstMatrix.add_result_vector(NaiveResultVector)" << endl;
        firstMatrix.create_comparison();
        cout << endl << "firstMatrix.create_comparison()" << endl;
        firstMatrix.show_the_matrix();
        cout << endl << "firstMatrix.show_the_matrix()" << endl;

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsedTime(end - start);
        cout << endl << elapsedTime.count() << " sec" << endl;

        return 0;
}
