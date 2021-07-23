
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



using namespace std;




int main(int argc, char** argv)
{
        string theFasta("10Bact.fa");
        //vector<string> allSequences(parse_this_fasta("10Bact.fa"));
        //cout << allSequences[0] << endl;
        TestTable refTable(10);
        NaiveIndex firstIndex(10,8);
        for(uint32 sequenceNumber = 0; sequenceNumber < 10; sequenceNumber++)
        //{
          //refTable.record_sequence(allSequences[sequenceNumber],sequenceNumber);
          //firstIndex.add_sketch(firstIndex.compute_sketch(allSequences[sequenceNumber], firstIndex.kmerSize));
        //}
        return 0;
}
