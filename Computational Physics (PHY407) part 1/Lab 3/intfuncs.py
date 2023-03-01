#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np


def trapz_int(a, b, N, f):
    """Trapezoid rule of integration
    IN:
    a, b: [floats] lower and upper bounds of integration
    N: [int] number of slices
    f: [function] a function to integrate
    OUT:
    [float] the integral
    PSEUDO-CODE:
    1. from a, b and N, compute h the int step
    2. initialize the result by adding (beginning + end)/2.
    3. loop to add the other terms f(a+kh)
    4. multiply the result by h to create final result
    """
    # 1.: from a, b and N, compute h
    h = (b-a)/N

    # 2.: the beginning and end
    result = 0.5*(f(a) + f(b))

    # 3. Then, loop over interior points
    for k in range(1, N):
        result += f(a + k*h)

    return h*result  # 4.


def simps_int(a, b, N, f):
    """Simpson's rule of integration
    IN:
    a, b: [floats] lower and upper bounds of integration
    N: [int] number of slices, HAS TO BE EVEN
    f: [function] a function to integrate
    OUT:
    [float] the integral
    PSEUDO-CODE:
    1. test if N is odd or even, and throw an error if odd
    2. from a, b and N, compute h the int step
    3. initialize the result by adding (beginning + end).
    4. loop over the odd k's to sum f(a+kh) from 1 to N-1
    5. add 4*(sum above) to the result
    6. loop over the even k's to sum f(a+kh) from 2 to N-2
    7. add 2*(sum above) to the result
    8. multiply the result by h/3 to create final result
    """

    # 1. test if N is odd or even, and throw an error if odd
    if N % 2 == 1:  # if True, the number is even
        # this aborts the calculation
        raise NameError("Simpson's rule needs an even number of slices!")

    # 2.: from a, b and N, compute h
    h = (b-a)/N

    # 3.: the beginning and end
    result = f(a) + f(b)

    # 4. loop over the odd k's to sum f(a+kh) from 1 to N-1
    odds = 0.
    for k in range(1, N, 2):
        odds += f(a + k*h)
    result += 4.*odds  # 5. add 4*(sum above) to the result

    # 6. loop over the even k's to sum f(a+kh) from 2 to N-2
    evens = 0.
    for k in range(2, N, 2):
        evens += f(a + k*h)
    result += 2*evens  # 7. add 2*(sum above) to the result

    return h*result/3  # 8. multiply the result by h/3 to create final result


def gauss_int(a, b, N, f):
    """ This function implements the Guassian quadrature integration from -1 to
    +1 in a similar way as the trapz and Simpson's rules above
    IN:
    a, b: [floats] lower and upper bounds of integration
    N: [int] number of sample points (NOT slices!)
    f: [function] a function to integrate
    new [bool]: True if we need to compute x, w, False otherwise
    OUT:
    [float] the integral"""
    x, w = gaussxwab(N, a, b)
    result = 0.
    for i in range(N):
        result += w[i]*f(x[i])
    return result


def gaussxw(N):
    """Newman's function for Gaussian quad from -1 to +1
    N is now the number of samplebegin{} points"""
    # Initial approximation to roots of the Legendre polynomial
    a = np.linspace(3, 4*N-1, N)/(4*N+2)
    x = np.cos(np.pi*a + 1/(8*N*N*np.tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta > epsilon:
        p0 = np.ones(N, float)
        p1 = np.copy(x)
        for k in range(1, N):
            p0, p1 = p1, ((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x, w


def gaussxwab(N, a, b):
    """Newman's function for Gaussian quad from a to b
    N is now the number of sample points"""
    x, w = gaussxw(N)
    return 0.5*(b-a)*x + 0.5*(b+a), 0.5*(b-a)*w

