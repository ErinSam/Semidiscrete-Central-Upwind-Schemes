"""
    The minmod function that is required by the van leer limiter
"""

# Erin Sam Joe | NITK '23 | Mechanical Engg. Major | Electronics & Electrical Engg. Minor

import numpy as np


def minmodScalar(var1, var2, var3):
    """
        Function that calculated the minmod of 3 variables
    """
    if ((var1<0) & (var2<0) & (var3<0)):
        return max(var1, var2, var3)
    elif ((var1>0) & (var2>0) & (var3>0)):
        return min(var1, var2, var3)
    else:
        return 0

def minmod(var1, var2, var3):
    """
        Function that calculated the minmod for 2 4-element vectors
    """
    mmd = np.zeros(4)

    for i in range(4):
        mmd[i] = minmodScalar(var1[i], var2[i], var3[i])

    return mmd

