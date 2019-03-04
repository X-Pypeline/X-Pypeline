# distutils: language = c++
# distutils: sources = xpipeline/cluster/src/fastsparseclusterprop.cpp

import numpy as np
cimport numpy as np
from libcpp cimport bool

cdef extern from "fastsparseclusterprop.h" nogil:
    void fastsparseclusterprop(np.double_t *labelledMap, np.double_t *likelihoodMap,
                               np.double_t *pixTime, np.double_t *pixFreq,
                               double clusterArray[],
                               const bool tf_properties, np.double_t *dimArray,
                               const int nClusters)

def clusterproperities_wrapper(object labelled_map, object likelihood,
                               tf_properties=False, pixel_time=None,
                               pixel_freq=None):
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
        np.ndarray[np.double_t, ndim=3, mode='c'] likelihood_tmp
        np.ndarray[np.double_t, ndim=1, mode='c'] pixel_time_tmp
        np.ndarray[np.double_t, ndim=1, mode='c'] pixel_freq_tmp
        np.ndarray[np.double_t, ndim=2, mode='c'] cluster_array

    # FIXME: This crazy reshaping is from the original logic being in MATLAB
    # I am sure with the numpy array format a cpp there is a way to not have to do this
    likelihood = likelihood.reshape(labelled_map.size, 1, -1)

    # FIXME the variable seems like it should not need to be passed!
    dimension_array_tmp = np.ascontiguousarray(likelihood.shape, dtype=np.float64)

    # If we do want tf_properities we need to have more columns in our array where
    # those values will go. If we do not then we just want the sum of likelihoods over clusters
    if tf_properties:
        cluster_array = np.zeros((7 + likelihood.shape[-1], labelled_map.max()), dtype=np.double)
    else:
        cluster_array = np.zeros((likelihood.shape[-1], labelled_map.max()), dtype=np.double)

    # Pass in the nearest neighbor labels of all clusters
    labelled_map_tmp = np.ascontiguousarray(labelled_map, dtype=np.float64)

    # Pass in a flattened array of likelihoods to calc statistics on
    likelihood_tmp = np.ascontiguousarray(likelihood, dtype=np.float64)

    # Pass in the start time of pixels
    if tf_properties:
        pixel_time_tmp = np.ascontiguousarray(pixel_time, dtype=np.float64)
    else:
        pixel_time_tmp = np.ascontiguousarray(np.zeros(labelled_map.size), dtype=np.float64)

    # pass in the end time of pixels
    if tf_properties:
        pixel_freq_tmp = np.ascontiguousarray(pixel_freq, dtype=np.float64)
    else:
        pixel_freq_tmp = np.ascontiguousarray(np.zeros(labelled_map.size), dtype=np.float64)

    fastsparseclusterprop(&labelled_map_tmp[0], &likelihood_tmp[0,0,0],
                          &pixel_time_tmp[0], &pixel_freq_tmp[0],
                          &cluster_array[0,0], tf_properties,
                          &dimension_array_tmp[0], labelled_map.max().astype(int))

    return cluster_array
