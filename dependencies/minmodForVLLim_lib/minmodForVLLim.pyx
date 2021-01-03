"""
    The minmod function that is required by the van leer limiter

    Sped up with cython
"""

# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor

import numpy as np

cimport numpy as np
cimport cython


cdef minmodScalar(double var1, double var2, double var3):
    """
        Function that calculated the minmod of 3 variables
    """
    if ((var1<0) & (var2<0) & (var3<0)):
        return max(var1, var2, var3)
    elif ((var1>0) & (var2>0) & (var3>0)):
        return min(var1, var2, var3)
    else:
        return 0


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef minmod(np.ndarray[double, ndim=1] var1, 
               np.ndarray[double, ndim=1] var2, 
               np.ndarray[double, ndim=1] var3):
    """
        Function that calculated the minmod for 2 4-element vectors
    """
    cdef np.ndarray[double, ndim=1] mmd = np.zeros(4)

    cdef int i
    for i in range(4):
        mmd[i] = minmodScalar(var1[i], var2[i], var3[i])

    return mmd


