# distutils: language = c++
# distutils: sources = xpipeline/cluster/src/fastsupercluster.cpp

import numpy as np
cimport numpy as np

cdef extern from "fastsupercluster.h" nogil:
     void fastsupercluster(const int N, np.double_t * rectangles, np.double_t * uncoveredMask)

def fastsupercluster_wrapper(ntriggers, object rectangles,):
    cdef:
        np.ndarray[np.double_t, ndim=2, mode='c'] rectangles_tmp
        np.ndarray[np.double_t, ndim=1, mode='c'] mask = np.zeros((ntriggers,), dtype=np.double)

    rectangles_tmp = np.ascontiguousarray(rectangles, dtype=np.float64)

    fastsupercluster(ntriggers, &rectangles_tmp[0,0], &mask[0])
    return mask
