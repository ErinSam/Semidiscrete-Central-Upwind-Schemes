"""
    Python script to convert from conserved variable vector to primitive variable vector and 
    vice versa

    Sped up with cython
"""

# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor


import numpy as np

cimport numpy as np
cimport cython


cpdef conservedToPrimitive(np.ndarray[double, ndim=1] U):
    """
        Function that converts a vector of conserved variables to vector of primitive variables.
        Meant for 2D compressible Euler equations. 

        Args:
            U (ndarray(4,)) : conserved variables vector

        Returns:
            W (ndarray(4,)) : primitive variables vector
    """
    cdef np.ndarray[double, ndim=1] W = np.zeros(4)
    cdef double c_ratio = 1.4
    
    W[0] = U[0]
    W[1] = U[1] / U[0]
    W[2] = U[2] / U[0] 
    W[3] = (c_ratio-1) * (U[3] - 0.5*U[0]*(U[1]**2 + U[2]**2))

    return W


cpdef primitiveToConserved(np.ndarray[double, ndim=1] W):
    """
        Function that converts a vector of primitive variables to a vector of conserved variables
        Meant for 2D compressible Euler equations. 

        Args:
            W (ndarray(4,)) : primitive variables vector

        Returns:
            U (ndarray(4,)) : conserved variables vector
    """
    cdef np.ndarray[double, ndim=1] U = np.zeros(4)
    cdef double c_ratio = 1.4

    U[0] = W[0]
    U[1] = W[0]*W[1]
    U[2] = W[0]*W[2] 
    U[3] = W[3]/(c_ratio-1) + 0.5*W[0]*(W[1]**2 + W[2]**2)

    return U
