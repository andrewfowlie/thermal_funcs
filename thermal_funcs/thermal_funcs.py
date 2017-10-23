"""
Wrapper for thermal functions from C++ library with SWIG
========================================================

There are various implmentations of the thermal functions J_F(y^2) and J_B(y^2), 
and their first- and second-derivatives.

They should be called via

>>> J_B(100., method='quad')
-0.0021509975248131265

>>> J_F(100., method='bessel')
0.002150965876834246

The method is an optional argument. First and second derivatives are
implemented by a keyword argument:

>>> J_F(100., derivative=1)
-9.324239650463552e-05

Tests:

>>> import numpy as np
>>> y_squared = np.linspace(-10, 10, 10)

>>> for method in METHODS:
...     np.vectorize(J_B)(y_squared, method)
...     np.vectorize(J_F)(y_squared, method)
array([-1.03301551, -2.37643967, -3.22920627, -3.4586863 , -2.87437545,
       -1.65904686, -1.13857283, -0.83787385, -0.63967438, -0.50041956])
array([ 2.55378622,  4.62466887,  4.57599995,  3.66746647,  2.44541129,
        1.53766395,  1.0938664 ,  0.81725016,  0.62904789,  0.49454272])
array([-1.03301518, -2.37643835, -3.22920637, -3.45867784, -2.87436571,
       -1.65904681, -1.13857282, -0.83787385, -0.63967438, -0.50041954])
array([ 2.55328922,  4.62461661,  4.57599754,  3.66745192,  2.44540827,
        1.53766392,  1.0938664 ,  0.81725016,  0.62904789,  0.4945427 ])
array([-9.7084274 , -7.62164949, -5.90432495, -4.42150614, -2.98134396,
       -1.7660153 , -2.10137679, -3.51278671, -5.88376321, -9.17182574])
array([ 15.43326144,  11.3026159 ,   7.6896565 ,   4.71314561,
         2.55512445,   1.64214217,   1.99504072,   3.22861748,
         5.20114362,   7.83237021])
array([ 8.2744023 ,  6.85296131,  5.32454196,  3.62990861,  1.59240946,
        1.59240946,  3.62990861,  5.32454196,  6.85296131,  8.2744023 ])
array([ 8.2744023 ,  6.85296131,  5.32454196,  3.62990861,  1.59240946,
        1.59240946,  3.62990861,  5.32454196,  6.85296131,  8.2744023 ])
array([-4.87947524, -5.29925691, -4.53529364, -2.66696586, -0.36007938,
       -0.47270639, -0.49809179, -0.42950048, -0.35893197, -0.29833225])
array([ 4.87947524,  5.29925691,  4.53529364,  2.66696586,  0.36007938,
        0.47270639,  0.49809179,  0.42950048,  0.35893197,  0.29833225])
array([-4.24227668, -4.52683182, -4.01645497, -2.66394847, -0.61658766,
       -0.47270639, -0.49809179, -0.42950048, -0.35893197, -0.29833225])
array([ 6.41771861,  6.83162332,  4.9203404 ,  2.34920927,  0.18597062,
        0.44674405,  0.48467656,  0.42254647,  0.35511487,  0.29613348])

Test derivative against numerical derivatives.

>>> np.vectorize(J_F)(y_squared, derivative=1, method="bessel") - np.vectorize(J_F)(y_squared, derivative=1, method="approx")
array([ -5.06825245e-03,   8.86800478e-05,   2.29852685e-05,
        -3.26269292e-04,  -2.71182328e-04,  -1.72205933e-07,
        -8.12912029e-08,   1.40641963e-07,  -4.34322898e-08,
         2.84767146e-08])

>>> np.vectorize(J_F)(y_squared, derivative=2, method="bessel") - np.vectorize(J_F)(y_squared, derivative=2, method="approx")
    array([ -5.28890689e-01,   4.34476713e-04,   2.28543869e-02,
             8.60189656e-03,   3.08677351e-02,   3.05994751e-02,
             1.99819211e-03,   4.25091227e-04,   1.31940097e-04,
             4.99716208e-05])

>>> np.vectorize(J_B)(y_squared, derivative=1, method="bessel") - np.vectorize(J_B)(y_squared, derivative=1, method="approx")
array([ -2.16951985e-04,  -2.55114955e-04,   6.97032234e-05,
        -3.38789396e-04,  -2.70543819e-04,   3.85995679e-07,
         1.32647813e-08,   1.30533110e-07,   1.23921645e-07,
        -3.70553466e-08])
        
>>> np.vectorize(J_B)(y_squared, derivative=2, method="bessel") - np.vectorize(J_B)(y_squared, derivative=2, method="approx")
array([ -1.18211875e-02,  -1.07143055e-02,  -2.78539577e-03,
         1.84933021e-02,   1.84354041e-01,   5.28891196e-02,
         2.65086492e-03,   5.05682480e-04,   1.48106517e-04,
         5.41656245e-05])
"""


# Load library


import _thermal_funcs


METHODS = ['quad', 'bessel', 'taylor', 'lim', 'approx', 'zeta']
METHODS_DERIV = ['bessel', 'approx']


def J_F(y_squared, method="bessel", derivative=0, **kwargs):
    """
    :param y_squared: Thermal function argument
    :type y_squared: float
    :param method: Method of calculating thermal function
    :type method: str
    :param derivative: Order of deriviatve
    :type derivative: int

    :returns: Thermal function J_F(y^2)
    :rtype: float
    """
    if derivative == 0:
        assert method in METHODS, "{} is not a known method: {}".format(method, METHODS)
        name = "J_F_{}".format(method)
    elif derivative == 1:
        assert method in METHODS_DERIV, "{} is not a known method: {}".format(method, METHODS_DERIV)
        name = "D1_J_F_{}".format(method)
    elif derivative == 2:
        assert method in METHODS_DERIV, "{} is not a known method: {}".format(method, METHODS_DERIV)
        name = "D2_J_F_{}".format(method)
    else:
        raise RuntimeError("Order must be 0 <= derivative <= 2")

    func = getattr(_thermal_funcs, name)
    return func(y_squared, **kwargs)


def J_B(y_squared, method="bessel", derivative=0, **kwargs):
    """
    :param y_squared: Thermal function argument
    :type y_squared: float
    :param method: Method of calculating thermal function
    :type method: str
    :param derivative: Order of deriviatve
    :type derivative: int

    :returns: Thermal function J_B(y^2)
    :rtype: float
    """
    if derivative == 0:
        assert method in METHODS, "{} is not a known method: {}".format(method, METHODS)
        name = "J_B_{}".format(method)
    elif derivative == 1:
        assert method in METHODS_DERIV, "{} is not a known method: {}".format(method, METHODS_DERIV)
        name = "D1_J_B_{}".format(method)
    elif derivative == 2:
        assert method in METHODS_DERIV, "{} is not a known method: {}".format(method, METHODS_DERIV)
        name = "D2_J_B_{}".format(method)
    else:
        raise RuntimeError("Order must be 0 <= derivative <= 2")

    func = getattr(_thermal_funcs, name)
    return func(y_squared, **kwargs)



if __name__ == '__main__':
    import doctest
    doctest.testmod(optionflags=doctest.REPORT_NDIFF)
    doctest.testmod()

