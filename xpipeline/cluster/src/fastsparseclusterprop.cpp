// replacement for fastclusterprop.cpp (fast version of regionprops.m) 
// it uses the list of pixels representation that 
// input : lablledMap, likelihoodMap, pixTime [optional], pixFreq [optional]
// output: clusterArray
// for details see clusterTFmapNew 

#include "fastsparseclusterprop.h"


inline double max(double a,double b) 
{
  if(a > b) {return a;}
  return b;
}

inline double min(double a,double b) 
{
  if(a < b) {return a;}
  return b;
}

void fastsparseclusterprop(const double *dimArray, const double *labelledMap, const double *likelihoodMap, const double *pixTime, const double *pixFreq, double *clusterArray, const bool doTFprops)
{
  //Declaration
  int nDims;

  // check input variable, if time/frequency scale information is
  // available do produce the cluster time/frequency information
  int nTFcols;

  if (doTFprops)
    { 
      nTFcols = 7;
    }
  else
    {
      nTFcols = 0;
    }

  // Number of dimesnion and size
  nDims=dimArray[2];
  if (3 != nDims && 2 != nDims)  {
    printf("%s\n","Error the number of dimension for likelihoodMap is not 2 or 3");
    return ;
  }
  int nLikelihoods;
  int colLen=dimArray[0];
  int rowLen=dimArray[1];

  if( 3== nDims) {
    nLikelihoods=dimArray[2]; }
  else {
    nLikelihoods=1;}
    
  if ( rowLen > 1 )
    {
      std::cout << "More than column, this function supports only the list of pixel format which should have a single column for the full TF map" << std::endl;
    }

  int nClusters=0;
  for(int j=0;j<colLen;j++)
    nClusters=(int)max(double(nClusters),labelledMap[j]);
  
  for(int j=0;j<colLen;j++){
    int label=int(labelledMap[j])-1;
    if( -1 == label) {continue;}
    // skip TF properties if not requested
    if (doTFprops)
      {
	// compute min t
	if(clusterArray[(0*nClusters)+label]>0){
	  clusterArray[(0*nClusters)+label] = 
	    min(clusterArray[(0*nClusters)+label], pixTime[j]-0.5);
	}
	else {clusterArray[(0*nClusters)+label] = pixTime[j]-0.5;}
	//compute mean t
	clusterArray[(1*nClusters)+label] = 
	  ( clusterArray[(1*nClusters)+label] * clusterArray[(nTFcols*nClusters)+label] +
	    pixTime[j] * likelihoodMap[0 + j]) /
	  ( clusterArray[(nTFcols*nClusters)+label] + likelihoodMap[0 + j]);
	// compute max t
	clusterArray[(2*nClusters)+label] = 
	  max(clusterArray[(2*nClusters)+label], pixTime[j]+0.5);
	// compute min f;
	if(clusterArray[(3*nClusters)+label]>0){
	  clusterArray[(3*nClusters)+label] = 
	    min(clusterArray[(3*nClusters)+label], pixFreq[j]-0.5);
	}
	else {clusterArray[(3*nClusters)+label] = pixFreq[j]-0.5;}
	//compute mean f
	clusterArray[(4*nClusters)+label] = 
      ( clusterArray[(4*nClusters)+label] * clusterArray[(nTFcols*nClusters)+label] +
	pixFreq[j] * likelihoodMap[0 + j]) /
	  ( clusterArray[(nTFcols*nClusters)+label] + likelihoodMap[0 + j]);
	//compute max f
	clusterArray[(5*nClusters)+label] = 
	  max(clusterArray[(5*nClusters)+label], pixFreq[j]+0.5);
	//compute Area
	clusterArray[(6*nClusters)+label] = clusterArray[(6*nClusters)+label]+1;
      }
    //always compute sum for each likelihood
    for(int k=0;k<nLikelihoods;k++) {
      clusterArray[((nTFcols+k)*nClusters)+label]=
	clusterArray[((nTFcols+k)*nClusters)+label]+
	likelihoodMap[k*colLen*rowLen + j];
    }
    
  }

  return;
}
    
