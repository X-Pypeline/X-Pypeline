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

inline map<int,int> create_map(vector<int> key)
{
    map<int,int> m;
    for (int i = 0; i < key.size(); ++i) {
        m[key[i]] = i+1;
    }
    return m;
}

vector<double> fastsparseclusterprop(const double *labelledMap, const double *likelihoodMap, const double *pixTime, const double *pixFreq, const bool doTFprops, const double *dimArray, const int nClusters)
{

    int colLen=dimArray[0];
    int rowLen=dimArray[1];
    int nLikelihoods=dimArray[2];
    int k = 0; //FIXME this should be a specific likelihood that you want to sort on
    double energy_of_cluster[nClusters];

    for(int j=0;j<colLen;j++){
        int label=int(labelledMap[j])-1;
        if( -1 == label) {continue;}
        // compute all enrgies of cluster for supplied likelihood
        // we use these values to only retain the top X percent
        // of clusters
        energy_of_cluster[label]=
            energy_of_cluster[label]+
            likelihoodMap[k*colLen*rowLen + j];
        }

    // We need to create a mask for labelledMap based on the
    // loudest X percent of clusters

    // Find the indices of the top 10 percent of
    // clusters
    int idx[nClusters];
    iota(idx, idx + nClusters, 0);
    sort(idx, idx + nClusters,
          [&](int lhs, int rhs) { return energy_of_cluster[lhs] < energy_of_cluster[rhs]; });
    int percentile_index = int(ceil(nClusters*0.1));
    vector<int> idx_of_loudest_clusters(idx, idx + percentile_index);
    sort(begin(idx_of_loudest_clusters), end(idx_of_loudest_clusters));

    // Map surviving clusters to "new" clusters starting at 1
    map<int,int> m = create_map(idx_of_loudest_clusters);

    // Prepare to reate a mask by looping over
    // These clusters and making a mask
    // of 0 if labelledMap[j] not in top 1 percent
    // of cluster and 1 if it is
    vector<int> mask(colLen, 0);

    for (int i = 0; i < colLen; ++i) {
        if (binary_search(begin(idx_of_loudest_clusters),
                          end(idx_of_loudest_clusters),
                          labelledMap[i])) {
            mask[i] = m[labelledMap[i]];
        }
    }

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

    vector<double> clusterArray((nTFcols + nLikelihoods)*percentile_index, 0);

    int counter = 0;
    if (doTFprops)
        {for(int j=0;j<colLen;j++){
            int label=int(mask[j])-1;
            if( -1 == label) {continue;}

            // compute min t
            if(clusterArray[(0*percentile_index)+label]>0){
              clusterArray[(0*percentile_index)+label] =
                min(clusterArray[(0*percentile_index)+label], pixTime[j]-0.5);
            }
            else {clusterArray[(0*percentile_index)+label] = pixTime[j]-0.5;}

            //compute mean t
            clusterArray[(1*percentile_index)+label] =
              ( clusterArray[(1*percentile_index)+label] * clusterArray[(nTFcols*percentile_index)+label] +
                pixTime[j] * likelihoodMap[0 + j]) /
              ( clusterArray[(nTFcols*percentile_index)+label] + likelihoodMap[0 + j]);

            // compute max t
            clusterArray[(2*percentile_index)+label] =
              max(clusterArray[(2*percentile_index)+label], pixTime[j]+0.5);

            // compute min f;
            if(clusterArray[(3*percentile_index)+label]>0){
              clusterArray[(3*percentile_index)+label] =
                min(clusterArray[(3*percentile_index)+label], pixFreq[j]-0.5);
            }
            else {clusterArray[(3*percentile_index)+label] = pixFreq[j]-0.5;}

            //compute mean f
            clusterArray[(4*percentile_index)+label] =
          ( clusterArray[(4*percentile_index)+label] * clusterArray[(nTFcols*percentile_index)+label] +
            pixFreq[j] * likelihoodMap[0 + j]) /
              ( clusterArray[(nTFcols*percentile_index)+label] + likelihoodMap[0 + j]);

            //compute max f
            clusterArray[(5*percentile_index)+label] =
              max(clusterArray[(5*percentile_index)+label], pixFreq[j]+0.5);

            //compute Area
            clusterArray[(6*percentile_index)+label] = clusterArray[(6*percentile_index)+label]+1;

            //always compute sum for each likelihood
            for(int k=0;k<nLikelihoods;k++) {
              clusterArray[((nTFcols+k)*percentile_index)+label]=
                clusterArray[((nTFcols+k)*percentile_index)+label]+
                likelihoodMap[k*colLen*rowLen + j];
            }
        }}
    else
        // Just return the likelihood values for the clusters
        {for(int j=0;j<colLen;j++){
            int label=int(mask[j])-1;
            if( -1 == label) {continue;}
            //cout << label <<endl;
            //always compute sum for each likelihood
            for(int k=0;k<nLikelihoods;k++) {
                clusterArray[((nTFcols+k)*percentile_index)+label]=
                    clusterArray[((nTFcols+k)*percentile_index)+label]+
                    likelihoodMap[k*colLen*rowLen + j];
            }
        }}

  // We need to create a mask for labelledMap based on the
  // loudest X percent of clusters

  // Find the indices of the top 10 percent of
  // clusters
  //int idx[nClusters];
  //iota(idx, idx + nClusters, 0);
  //sort(idx, idx + nClusters,
  //        [&](int lhs, int rhs) { return clusterArray[lhs] < clusterArray[rhs]; });
  //int percentile_index = int(ceil(nClusters*0.1));
  //unordered_set<int> idx_of_loudest_clusters(idx, idx + percentile_index);

  // Prepare to reate a mask by looping over
  // These clusters and making a mask
  // of 0 if labelledMap[j] not in top 1 percent
  // of cluster and 1 if it is
  //int mask[colLen]

  //for (int i = 0; i < colLen; ++i) {
  //    if (idx_of_loudest_clusters.find(labelledMap[i]) != idx_of_loudest_clusters.end()) {
  //        mask[i] = 1;
  //    else
  //        mask[i] = 0;
  //    }
  //}
  return clusterArray;
}
