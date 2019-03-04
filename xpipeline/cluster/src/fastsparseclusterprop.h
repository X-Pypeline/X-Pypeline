#include "math.h"
#include "stdio.h"
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

void fastsparseclusterprop(const double *labelledMap, const double * likelihoodMap, const double * pixTime, const double * pixFreq, double clusterArray[], const bool doTFprops, const double *dimArray, const int nClusters);
