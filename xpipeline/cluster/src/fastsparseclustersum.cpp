// replacement for fastclustersum.cpp (fast version of regionsums.m) 
// it uses the list of pixels representation that 
// input : lablledMap, likelihoodMap, pixTime [optional], pixFreq [optional]
// output: clusterArray
// for details see clusterTFmapNew 

#include "fastsparseclustersum.h"


void fastsparseclustersum(const int colLen, const int nClusters, const double *labelledMap, const double *likelihoodMap, double *clusterArray)
{
  //Declaration
  int nTFcols;
  int nLikelihoods;
  int rowLen;
  nLikelihoods = 1;
  rowLen = 1;
  nTFcols = 0;

  for(int j=0;j<colLen;j++){
    int label=int(labelledMap[j])-1;
    if( -1 == label) {continue;}
    //compute sum for each likelihood
    for(int k=0;k<nLikelihoods;k++) {
      clusterArray[((nTFcols+k)*nClusters)+label]=
	clusterArray[((nTFcols+k)*nClusters)+label]+
	likelihoodMap[k*colLen*rowLen + j];
    }
    
  }

  return;
}
