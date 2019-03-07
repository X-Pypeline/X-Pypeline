# distutils: language = c++
# distutils: sources = xpipeline/cluster/src/fastsparseclusterprop.cpp

import numpy
cimport numpy
from libcpp cimport bool
from libcpp.vector cimport vector

cdef extern from "fastsparseclusterprop.h" nogil:
    vector[double] fastsparseclusterprop(numpy.double_t *labelledMap, numpy.double_t *likelihoodMap,
                          numpy.double_t *pixTime, numpy.double_t *pixFreq,
                          const bool tf_properties, numpy.double_t *dimArray,
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
        numpy.ndarray[numpy.double_t, ndim=1, mode='c'] dimension_array_tmp
        numpy.ndarray[numpy.double_t, ndim=1, mode='c'] labelled_map_tmp
        numpy.ndarray[numpy.double_t, ndim=3, mode='c'] likelihood_tmp
        numpy.ndarray[numpy.double_t, ndim=1, mode='c'] pixel_time_tmp
        numpy.ndarray[numpy.double_t, ndim=1, mode='c'] pixel_freq_tmp

    # FIXME: This crazy reshaping is from the original logic being in MATLAB
    # I am sure with the numpy array format a cpp there is a way to not have to do this
    likelihood = likelihood.reshape(labelled_map.size, 1, -1)

    # FIXME the variable seems like it should not need to be passed!
    dimension_array_tmp = numpy.ascontiguousarray(likelihood.shape, dtype=numpy.float64)

    # Pass in the nearest neighbor labels of all clusters
    labelled_map_tmp = numpy.ascontiguousarray(labelled_map, dtype=numpy.float64)

    # Pass in a flattened array of likelihoods to calc statistics on
    likelihood_tmp = numpy.ascontiguousarray(likelihood, dtype=numpy.float64)

    # Pass in the start time and start frequency of pixels
    if tf_properties:
        pixel_time_tmp = numpy.ascontiguousarray(pixel_time, dtype=numpy.float64)
        pixel_freq_tmp = numpy.ascontiguousarray(pixel_freq, dtype=numpy.float64)
        number_of_columns = int(dimension_array_tmp[2] + 7)
    else:
        pixel_freq_tmp = numpy.ascontiguousarray(numpy.zeros(labelled_map.size), dtype=numpy.float64)
        pixel_time_tmp = numpy.ascontiguousarray(numpy.zeros(labelled_map.size), dtype=numpy.float64)
        number_of_columns = int(dimension_array_tmp[2])

    cluster_array = fastsparseclusterprop(&labelled_map_tmp[0], &likelihood_tmp[0,0,0],
                          &pixel_time_tmp[0], &pixel_freq_tmp[0],
                          tf_properties, &dimension_array_tmp[0], labelled_map.max().astype(int))

    return numpy.asarray(cluster_array).reshape(number_of_columns, -1).T 
