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
array([-1.03301539, -2.37643826, -3.22920637, -3.45868482, -2.87437645,
       -1.65904685, -1.13857283, -0.83787385, -0.63967438, -0.50041956])
array([ 2.55328859,  4.62471422,  4.5760003 ,  3.66746892,  2.44541164,
        1.53766395,  1.09386639,  0.81725016,  0.62904789,  0.49454272])
array([-9.70842739, -7.62164956, -5.90432496, -4.42150614, -2.98134396,
       -1.7660153 , -2.10137678, -3.51278666, -5.88376315, -9.1718257 ])
array([ 15.43326181,  11.30261604,   7.68965655,   4.71314562,
         2.55512445,   1.64214217,   1.99504074,   3.22861748,
         5.20114365,   7.83237033])
array([ 9.45467251,  7.83047555,  6.08404071,  4.14768293,  1.81955257,
        1.81955257,  4.14768293,  6.08404071,  7.83047555,  9.45467251])
array([ 9.45467251,  7.83047555,  6.08404071,  4.14768293,  1.81955257,
        1.81955257,  4.14768293,  6.08404071,  7.83047555,  9.45467251])
array([-4.87947524, -5.29925691, -4.53529364, -2.66696586, -0.36007938,
       -0.47270639, -0.49809179, -0.42950048, -0.35893197, -0.29833225])
array([ 4.87947524,  5.29925691,  4.53529364,  2.66696586,  0.36007938,
        0.47270639,  0.49809179,  0.42950048,  0.35893197,  0.29833225])
array([-4.24227668, -4.52683182, -4.01645497, -2.66394847, -0.61658766,
       -0.47270639, -0.49809179, -0.42950048, -0.35893197, -0.29833225])
array([ 6.41771861,  6.83162332,  4.9203404 ,  2.34920927,  0.18597062,
        0.44674405,  0.48467656,  0.42254647,  0.35511487,  0.29613348])
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

