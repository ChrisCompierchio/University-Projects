# SolveLinear.py
# Python module for PHY407
# Paul Kushner, 2015-09-26
# Modifications by Nicolas Grisouard, 2018-09-26
# This module contains useful routines for solving linear systems of equations.
# Based on gausselim.py from Newman
from numpy import empty
import numpy as np
from numpy.random import rand
# The following will be useful for partial pivoting
from numpy import empty, copy

def PartialPivot(A_in, v_in):
    for i in range(len(v_in)-1):
        n = A_in[i]
        if np.abs(n[i]) < np.abs(n[i+1]):
            n = A_in[i+1]
            A_in[0], A_in[i] = copy(A_in[i]), copy(A_in[0])
            v_in[0], v_in[i] = copy(v_in[i]), copy(v_in[0]) 
    return A_in, v_in

def GaussElim(A_in, v_in):
    """Implement Gaussian Elimination. This should be non-destructive for input
    arrays, so we will copy A and v to
    temporary variables
    IN:
    A_in, the matrix to pivot and triangularize
    v_in, the RHS vector
    OUT:
    x, the vector solution of A_in x = v_in """
    # copy A and v to temporary variables using copy command
    A = copy(A_in)
    v = copy(v_in)
    N = len(v)

    for m in range(N):
        # Divide by the diagonal element
        div = A[m, m]
        if div == 0:
            PartialPivot(A, v)
        else:
            A[m, :] /= div
            v[m] /= div

        # Now subtract from the lower rows
        for i in range(m+1, N):
            mult = A[i, m]
            A[i, :] -= mult*A[m, :]
            v[i] -= mult*v[m]

    # Backsubstitution
    # create an array of the same type as the input array
    x = empty(N, float)
    for m in range(N-1, -1, -1):
        x[m] = v[m]
        for i in range(m+1, N):
            x[m] -= A[m, i]*x[i]
    return x


