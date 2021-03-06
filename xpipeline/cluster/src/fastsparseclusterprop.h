#include "math.h"
#include "stdio.h"
#include <string>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <vector>
#include <map>
#include <cmath>
#include <iterator>

using namespace std;

vector<double> fastsparseclusterprop(const double *labelledMap, const double * likelihoodMap, const double * pixTime, const double * pixFreq, const bool doTFprops, const double *dimArray, const int nClusters, const double *projectedAsdMagnitudeSquared);
