
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

        string theFasta("fastaest.fa");
        cout << endl << "the fasta" << endl;
        TestTable refTable(6);//genomes number
        cout << endl << "reftable" << endl;
        NaiveIndex firstIndex(512,6,8);//((bucketnumbers, genomes number ,bits kept for record))
        cout << endl << "firstindex" << endl;
        refTable.parse_fasta_for_refTable(theFasta);
        cout << endl << "parse_fasta_for_refTable" << endl;
        firstIndex.index_sequences_from_fasta(theFasta);
        cout << endl << "add_fasta_for_naive" << endl;
        string otherFasta("fastest.fa");
        cout << endl << "otherFasta" << endl;
        ifstream oneRead(otherFasta);
        cout << endl << "oneRead" << endl;
        string oneGenome("atgggcgggcgggcggggctttattagctga");
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

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsedTime(end - start);
        cout << endl << elapsedTime.count() << " sec" << endl;

        return 0;
}
