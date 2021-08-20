
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
#include <math.h>
#include <cassert>
#include "NaiveIndex.h"
#include "EvaluationMaker.h"
#include "TestTable.h"






int main(int argc, char** argv)
{
        //time
        std::chrono::time_point<std::chrono::system_clock> start, end, startDB, endDB, startQuery, endQuery;
        start = std::chrono::system_clock::now();

        //create ComparisonMatrix
        ComparisonMatrix firstMatrix;


        //production of the (hash) Testtable
        string theFasta("extract_from_10B.fa");
        TestTable refTable;
        refTable.parse_fasta_for_refTable(theFasta);
        ifstream theFileQuery("little_extract_from_extract_from10B.fa");
        string oneGenome(refTable.get_line_fasta_for_testtable(&theFileQuery));
        vector<long double> TestResultVector(refTable.query_belonging_genome(oneGenome));
        vector<pair<long double,uint16>> sortedVectorTest(refTable.sort_scores(TestResultVector,0));
        refTable.show_sorted_scores(sortedVectorTest);
        firstMatrix.add_result_vector(TestResultVector);

        //production of the index(es)
        for (uint binaryPower = 9; binaryPower < 25; binaryPower++)
        {
        string theFasta("extract_from_10B.fa");
        assert(pow(2,binaryPower) > pow(2,8) && (pow(2,binaryPower) < pow(2,30)));
        startDB = std::chrono::system_clock::now();
        NaiveIndex firstIndex(pow(2,binaryPower),256);//((bucketnumbers, binary decimal (2,4,8...1024...) for record))
        firstIndex.index_sequences_from_fasta(theFasta);
        endDB = std::chrono::system_clock::now();
        ifstream theFileQuery("little_extract_from_extract_from10B.fa");
        startQuery = std::chrono::system_clock::now();
        string oneGenome(firstIndex.get_line_fasta_for_naive(&theFileQuery));
        vector<long double> NaiveResultVector(firstIndex.query_sequence(oneGenome));
        endQuery = std::chrono::system_clock::now();
        vector<pair<long double,uint16>> sortedVector(firstIndex.sort_scores(NaiveResultVector,0));
        firstIndex.show_sorted_scores(sortedVector);
        firstMatrix.add_result_vector(NaiveResultVector);
        firstMatrix.fill_time((binaryPower - 8), endDB-startDB, endQuery-startQuery);
        }

        //The comparison
        firstMatrix.create_comparison();
        firstMatrix.show_the_matrix();

        //the record chosen
        firstMatrix.write_result("Jaccard_index_result_with_time.csv", "jaccard", 1);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsedTime(end - start);
        cout << endl << elapsedTime.count() << " sec" << endl;

        return 0;
}
