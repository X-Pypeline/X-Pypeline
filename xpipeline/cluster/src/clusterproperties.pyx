# distutils: language = c++
# distutils: sources = xpipeline/cluster/src/fastsparseclusterprop.cpp

import numpy as np
cimport numpy as np

cdef extern from "fastsparseclusterprop.h" nogil:
    void fastsparseclusterprop(np.int_t *dimArray, np.double_t *labelledMap,
                               np.double_t * likelihoodMap, np.double_t * pixTime, 
                               np.double_t * pixFreq, np.double_t *clusterArray)

def clusterproperities_wrapper(object dimension_array, object labelled_map,
                               object likelihood_map,object pixel_time,
                               object pixel_freq):
    cdef:
        np.ndarray[np.int_t, ndim=1, mode='c'] dimension_array_tmp
        np.ndarray[np.double_t, ndim=1, mode='c'] labelled_map_tmp
        np.ndarray[np.double_t, ndim=1, mode='c'] coord_dim_array_tmp
        np.ndarray[np.double_t, ndim=1, mode='c'] pixel_time_tmp
        np.ndarray[np.double_t, ndim=1, mode='c'] pixel_freq_array_tmp
        np.ndarray[np.double_t, ndim=1, mode='c'] cluster_array = np.zeros((labelled_map.max, 7 + dimension_array[2]), dtype=np.double)

    labelled_map_tmp = np.ascontiguousarray(labelled_map, dtype=np.float64)
    likelihood_map_tmp = np.ascontiguousarray(likelihood_map, dtype=np.float64)
    pixel_time_tmp = np.ascontiguousarray(pixel_time, dtype=np.float64)
    pixel_freq_tmp = np.ascontiguousarray(pixel_freq, dtype=np.float64)
    

    fastsparseclusterprop(&dimension_array_tmp[0], &labelled_map_tmp[0],
                          &likelihood_map_tmp[0], &pixel_time_tmp[0], &pixel_freq_tmp[0])
    return cluster_array
