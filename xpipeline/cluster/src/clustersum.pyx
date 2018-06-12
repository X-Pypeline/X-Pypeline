# distutils: language = c++
# distutils: sources = xpipeline/cluster/src/fastsparseclustersum.cpp

import numpy as np
cimport numpy as np

cdef extern from "fastsparseclustersum.h" nogil:
    void fastsparseclustersum(const int colLen,
                              const int nClusters,
                              np.double_t *labelledMap,
                              np.double_t *likelihoodMap,
                              np.double_t *clusterArray)

def clustersum_wrapper(object labelled_map,
                       object likelihood_map):
    """This outputs properities of every cluster on the TF map
        Returns:
            cluster_array (array):
                mx8+ matrix (m is # of clusters)
                column 0: minimum time of cluster
                column 1: weighted center time of cluster
                column 2: maximum time of cluster
                column 3: minimum frequency of cluster
                column 4: weighted center frequency of cluster
                column 5: maximum frequency of cluster
                column 6: number of pixels in cluster
                column 7-?: sum-over-cluster map values for each likelihood
                All times and frequencies are in units of bins, unless
                mapDim specified.  For example, a single-pixel cluster
                at (time,freq) bin = (100,20) would have clusterArray
                values
                 [ 99.5 100 100.5 19.5 20 20.5 1 ... likelihoods ... ]
                If mapDim is specified, the output is in seconds and
                Hz.
    """
    cdef:
        np.ndarray[np.double_t, ndim=1, mode='c'] labelled_map_tmp
        np.ndarray[np.double_t, ndim=1, mode='c'] likelihood_map_tmp
        np.ndarray[np.double_t, ndim=1, mode='c'] cluster_array = np.zeros((labelled_map.max(),), dtype=np.double)

    labelled_map_tmp = np.ascontiguousarray(labelled_map, dtype=np.float64)
    likelihood_map_tmp = np.ascontiguousarray(likelihood_map, dtype=np.float64)
    col_len = likelihood_map_tmp.size
    n_clusters = labelled_map_tmp.max().astype(int)


    fastsparseclustersum(col_len, n_clusters,
                         &labelled_map_tmp[0],
                         &likelihood_map_tmp[0],
                         &cluster_array[0])
    return cluster_array
