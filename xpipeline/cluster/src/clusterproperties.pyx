# distutils: language = c++
# distutils: sources = xpipeline/cluster/src/fastsparseclusterprop.cpp

import numpy as np
cimport numpy as np

cdef extern from "fastsparseclusterprop.h" nogil:
    void fastsparseclusterprop(np.double_t *dimArray, np.double_t *labelledMap,
                               np.double_t *likelihoodMap, np.double_t *pixTime,
                               np.double_t *pixFreq, np.double_t *clusterArray)

def clusterproperities_wrapper(object dimension_array, object labelled_map,
                               object likelihood_map, object pixel_time,
                               object pixel_freq):
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
        np.ndarray[np.double_t, ndim=1, mode='c'] dimension_array_tmp
        np.ndarray[np.double_t, ndim=1, mode='c'] labelled_map_tmp
        np.ndarray[np.double_t, ndim=1, mode='c'] likelihood_map_tmp
        np.ndarray[np.double_t, ndim=1, mode='c'] pixel_time_tmp
        np.ndarray[np.double_t, ndim=1, mode='c'] pixel_freq_tmp
        np.ndarray[np.double_t, ndim=2, mode='c'] cluster_array = np.zeros((8, labelled_map.max()), dtype=np.double)

    dimension_array_tmp = np.ascontiguousarray(dimension_array, dtype=np.float64)
    labelled_map_tmp = np.ascontiguousarray(labelled_map, dtype=np.float64)
    likelihood_map_tmp = np.ascontiguousarray(likelihood_map, dtype=np.float64)
    pixel_time_tmp = np.ascontiguousarray(pixel_time, dtype=np.float64)
    pixel_freq_tmp = np.ascontiguousarray(pixel_freq, dtype=np.float64)


    fastsparseclusterprop(&dimension_array_tmp[0], &labelled_map_tmp[0],
                          &likelihood_map_tmp[0],
                          &pixel_time_tmp[0], &pixel_freq_tmp[0],
                          &cluster_array[0,0])
    return cluster_array
