"""
    CPython script for the defining functions to calculate the interface speed
"""

import numpy as np 
import sys

cimport numpy as np
cimport cython

sys.path.append("./../fluxRelated_lib")
from fluxRelated import eigvalMaxMinF, eigvalMaxMinG


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef x(np.ndarray[double, ndim=1] fieldEast, np.ndarray[double, ndim=1] fieldWest, 
        bint minim=False):
    """
        Function that calculates the interface speed of the right interface of the cell.

        Args:
            fieldEast (ndarray(4,)) : the east flowField of the owner cell
            fieldWest (ndarray(4,)) : the west flowField of the neighbor cell
            minim (bool) : the 

        Returns:
            intfSpeed (double) : the required interface speed
    """
    cdef double eigval1, eigval2, intfSpeed

    if ( minim ):
        # Case that we want the minimum interface speed
        eigval1 = eigvalMaxMinF(fieldWest, minim=True)
        eigval2 = eigvalMaxMinF(fieldEast, minim=True)
        intfSpeed = min(eigval1, eigval2, 0)
        return intfSpeed

    else: 
        # Case that we want the maximum interface speed
        eigval1 = eigvalMaxMinF(fieldWest)
        eigval2 = eigvalMaxMinF(fieldEast)
        intfSpeed = max(eigval1, eigval2, 0)
        return intfSpeed

        
@cython.boundscheck(False)
@cython.wraparound(False)
cpdef y(np.ndarray[double, ndim=1] fieldNorth, np.ndarray[double, ndim=1] fieldSouth, 
        bint minim=False):
    """
        Function that calculates the interface speed of the top interface of the cell.

        Args:
            fieldNorth (ndarray(4,)) : the north flowField of the owner cell
            fieldSouth (ndarray(4,)) : the south flowField of the neighbor cell
            minim (bool) : the 

        Returns:
            intfSpeed (double) : the required interface speed
    """
    cdef double eigval1, eigval2, intfSpeed

    if ( minim ):
        # Case that we want the minimum interface speed
        eigval1 = eigvalMaxMinG(fieldSouth, minim=True)
        eigval2 = eigvalMaxMinG(fieldNorth, minim=True)
        intfSpeed = min(eigval1, eigval2, 0) 
        return intfSpeed

    else: 
        # Case that we want the maximum interface speed
        eigval1 = eigvalMaxMinG(fieldSouth)
        eigval2 = eigvalMaxMinG(fieldNorth)
        intfSpeed = max(eigval1, eigval2, 0)
        return intfSpeed


