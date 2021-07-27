
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

        return 0;
}
