"""
Wrapper for thermal functions from C++ library with SWIG
========================================================

There are various implmentations of the thermal functions J_F(y^2) and J_B(y^2).

They should be called via

>>> J_B(100., method='quad')
-0.0021509975248131265

>>> J_F(100., method='bessel')
0.002150965877087245

The method is an optional argument.

>>> import numpy as np
>>> y_squared = np.linspace(-10, 10, 10)
>>> for method in METHODS:
...     np.vectorize(J_B)(y_squared, method)
...     np.vectorize(J_F)(y_squared, method)
array([-1.03301551, -2.37643967, -3.22920627, -3.4586863 , -2.87437545,
       -1.65904686, -1.13857283, -0.83787385, -0.63967438, -0.50041956])
array([ 2.55378622,  4.62466887,  4.57599995,  3.66746647,  2.44541129,
        1.53766395,  1.0938664 ,  0.81725016,  0.62904789,  0.49454272])
array([-1.03301551, -2.37643979, -3.22920623, -3.45868664, -2.87437645,
       -1.65904686, -1.13857283, -0.83787385, -0.63967438, -0.50041956])
array([ 2.55392427,  4.62465458,  4.57600028,  3.66746892,  2.44541164,
        1.53766395,  1.09386639,  0.81725016,  0.62904789,  0.49454272])
array([-9.69735519, -7.61783034, -5.90338519, -4.42139071, -2.98134261,
       -1.76601401, -2.10127724, -3.51205276, -5.88106505, -9.16475924])
array([ 16.05591746,  11.59697008,   7.79738032,   4.73650148,
         2.5559924 ,   1.64127158,   1.97147   ,   3.11922085,
         4.90027541,   7.19158449])
array([ 9.45467251,  7.83047555,  6.08404071,  4.14768293,  1.81955257,
        1.81955257,  4.14768293,  6.08404071,  7.83047555,  9.45467251])
array([ 9.45467251,  7.83047555,  6.08404071,  4.14768293,  1.81955257,
        1.81955257,  4.14768293,  6.08404071,  7.83047555,  9.45467251])
array([-4.87947524, -5.29925691, -4.53529364, -2.66696586, -0.36007938,
       -0.47270639, -0.49809179, -0.42950048, -0.35893197, -0.29833225])
array([ 4.87947524,  5.29925691,  4.53529364,  2.66696586,  0.36007938,
        0.47270639,  0.49809179,  0.42950048,  0.35893197,  0.29833225])
array([-4.24227668, -4.52683182, -4.01645497, -2.66394847, -0.61658766,
       -0.50629843, -0.51317766, -0.43694993, -0.36292338, -0.30059963])
array([ 6.41771861,  6.83162332,  4.9203404 ,  2.34920927,  0.18597062,
        0.44674406,  0.48467656,  0.42254647,  0.35511487,  0.29613349])
"""


# Load library


import _thermal_funcs


METHODS = ['quad', 'bessel', 'taylor', 'lim', 'approx', 'zeta']


def J_F(y_squared, method="bessel", **kwargs):
    """
    :param y_squared: Thermal function argument
    :type y_squared: float
    :param method: Method of calculating thermal function
    :type method: str

    :returns: Thermal function J_F(y^2)
    :rtype: float
    """
    try:
        name = "J_F_{}".format(method)
        func = getattr(_thermal_funcs, name)
    except AttributeError:
        error = "{} is not a known method: {}".format(method, METHODS)
        raise RuntimeError(error)

    return func(y_squared, **kwargs)


def J_B(y_squared, method="bessel", **kwargs):
    """
    :param y_squared: Thermal function argument
    :type y_squared: float
    :param method: Method of calculating thermal function
    :type method: str

    :returns: Thermal function J_B(y^2)
    :rtype: float
    """
    try:
        name = "J_B_{}".format(method)
        func = getattr(_thermal_funcs, name)
    except AttributeError:
        error = "{} is not a known method: {}".format(method, METHODS)
        raise RuntimeError(error)

    return func(y_squared, **kwargs)


if __name__ == '__main__':
    import doctest
    doctest.testmod()

