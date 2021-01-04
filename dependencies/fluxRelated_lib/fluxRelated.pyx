"""
    Python script that contains definitions of the functions that are required for 
    calculating the flux related values that do not require mesh topology information

    Sped up with cython
"""

# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor


import numpy as np
import math as m 

cimport numpy as np
cimport cython


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef eigvalMaxMinF(np.ndarray[double, ndim=1] U, bint minim=False):
    """
        Function that finds the maximum or the minimum eigenvalue (based on choice) of 
        the Jacobian matrix of the x-direction flux term for Compressible 2D Euler 
        Equations. Eigenvalues of \partial F/ \partial U. 

        Args:
            U (ndarray(4,)) : conserved flow field variables 
            minmax (str) : whether to provide max. or min. eigenvalue of the Jacobian

        Returns: 
            eval (float) : max. or min. eigenvalue of the Jacobian
    """
    # Heat capacity ratio
    cdef double c_ratio = 1.4
    cdef double c_ = 1-c_ratio

    # For my ease
    cdef double u1 = U[0]
    cdef double u2 = U[1]
    cdef double u3 = U[2]
    cdef double u4 = U[3]

    # Creating the Jacobian matrix
    cdef np.ndarray[double, ndim=2] jacb = np.array([[0, 1, 0, 0],
                     [(-u2**2/u1**2 + c_/2*(u2**2/u1**2 + u3**2/u1**2)), (2*u2/u1 - c_*u2/u1), (-c_*u3/u1), c_],
                     [(-u2*u3/u1**2), u3/u1, u2/u1, 0],
                     [(-c_ratio*u4*u2/u1**2 + c_*(u2**3/u1**3 + u3**3/u1**3)), (c_ratio*u4/u1 - 1.5*c_*u2**2/u1**2), (-1.5*c_*u3**2/u1**2), c_ratio*u2/u1]])

    # Obtaining the max or the min of the eigenvalue
    cdef np.ndarray[double complex, ndim=1] eigvals = np.linalg.eig(jacb)[0]
    if ( minim ):
        return eigvals.min().real
    else:
        return eigvals.max().real


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef eigvalMaxMinG(np.ndarray[double, ndim=1] U, bint minim=False):
    """
        Function that finds the maximum or the minimum eigenvalue (based on choice) of 
        the Jacobian matrix of the y-direction flux term for Compressible 2D Euler 
        Equations. Eigenvalues of \partial G/ \partial U. 

        Args:
            U (ndarray(4,)) : conserved flow field variables 
            minmax (str) : whether to provide max. or min. eigenvalue of the Jacobian

        Returns: 
            eval (float) : max. or min. eigenvalue of the Jacobian
    """

    # Heat capacity ratio
    cdef double c_ratio = 1.4
    cdef double c_ = 1-c_ratio

    # For my ease
    cdef double u1 = U[0]
    cdef double u2 = U[1]
    cdef double u3 = U[2]
    cdef double u4 = U[3]

    # Creating the Jacobian matrix
    cdef np.ndarray[double, ndim=2] jacb = np.array([[0, 0, 1, 0],
                                 [(-u2*u3/u1**2), u3/u1, u2/u1, 0],
                                 [(-u3**2/u1**2 + c_/2*(u2**2/u1**2 + u3**2/u1**2)), -c_*u3/u1, \
                                                                    (2*u3/u1 - c_*u3/u1), c_],
                                 [(-c_ratio*u4*u3/u1**2 + c_*(u2**3/u1**3 + u3**3/u1**3)), \
                                    (-1.5*c_*u2**2/u1**2), (c_ratio*u4/u1 - 1.5*c_*u3**2/u1**2), \
                                                                        c_ratio*u3/u1]])

    # Obtaining the max or the min of the eigenvalue
    cdef np.ndarray[double complex, ndim=1] eigvals = np.linalg.eig(jacb)[0]
    if ( minim ):
        return eigvals.min().real
    else:
        return eigvals.max().real


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef fluxF(np.ndarray[double, ndim=1] U):
    """
        Function that calculated the flux from a given conserved flowfield
        variable vector in the x-direction for Compressible 2D Euler Equations.
        F(U) from 
            U_t + F(U)_x + G(U)_y = 0

        Args:
            U (ndarray(4,)) : conserved flow field variables 

        Returns:
            flux (ndarray(4,)) : flux in the x-direction
    """
    # Heat Capacity Ratio
    cdef double c_ratio = 1.4

    # Creating the flux vector
    cdef np.ndarray[double, ndim=1] flux = np.zeros(4)

    # Calculating flux 
    flux[0] = U[1]
    flux[1] = pow(U[1],2)/U[0] + (c_ratio-1) * (U[3] - 0.5*pow(U[1],2)/U[0] - 0.5*pow(U[2],2)/U[0])
    flux[2] = U[1]*U[2]/U[0]
    flux[3] = U[1]/U[0] * (c_ratio*U[3] - (c_ratio-1) * (0.5*pow(U[1],2)/U[0] - \
                                                            0.5*pow(U[2],2)/U[0]))

    return flux


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef fluxG(np.ndarray[double, ndim=1] U):
    """
        Function that calculated the flux from a given conserved flowfield
        variable vector in the x-direction for Compressible 2D Euler Equations.
        G(U) from 
            U_t + F(U)_x + G(U)_y = 0

        Args:
            U (ndarray(4,)) : conserved flow field variables 

        Returns:
            flux (ndarray(4,)) : flux in the x-direction
    """
    # Heat Capacity Ratio
    cdef double c_ratio = 1.4

    # Creating the flux vector
    cdef np.ndarray[double, ndim=1] flux = np.zeros(4)

    # Calculating flux 
    flux[0] = U[2]
    flux[1] = U[1]*U[2]/U[0] 
    flux[2] = pow(U[2],2)/U[0] + (c_ratio-1) * (U[3] - 0.5*pow(U[1],2)/U[0] - 0.5*pow(U[2],2)/U[0]) 
    flux[3] = U[2]/U[0] * (c_ratio*U[3] - (c_ratio-1) * (0.5*pow(U[1],2)/U[0] - \
                                                            0.5*pow(U[2],2)/U[0]))

    return flux

