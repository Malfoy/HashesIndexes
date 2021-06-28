#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <atomic>
#include <mutex>
#include <stdint.h>
#include "NaiveIndex.h"


using namespace std;



NaiveIndex::NaiveIndex(int64_t n){
  matrix.resize(n);
  nb_minimizer=n;
  nb_genomes=0;
};
