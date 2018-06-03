# distutils: language = c++
# distutils: sources = xpipeline/cluster/src/fastlabel.cpp

import numpy as np
cimport numpy as np

cdef extern from "fastlabel.h" nogil:
     void fastlabel(const int nPixels, np.double_t * coords,
                    np.double_t * coordDim, const int nNeighboors,
                    np.double_t * labelList)

def fastlabel_wrapper(object coord_array, object coord_dim_array,
                      connectivity, npixels):
    cdef:
        np.ndarray[np.double_t, ndim=2, mode='c'] coord_array_tmp
        np.ndarray[np.double_t, ndim=1, mode='c'] coord_dim_array_tmp
        np.ndarray[np.double_t, ndim=1, mode='c'] label_list = np.zeros((npixels,), dtype=np.double)

    coord_array_tmp = np.ascontiguousarray(coord_array, dtype=np.float64)
    coord_dim_array_tmp = np.ascontiguousarray(coord_dim_array, dtype=np.float64)

    fastlabel(npixels, &coord_array_tmp[0,0], &coord_dim_array_tmp[0],
              connectivity, &label_list[0])
    return label_list
