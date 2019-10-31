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
>>> np.set_printoptions(precision=6)
>>> y_squared = np.linspace(-10, 10, 10)

Test thermal functions.

>>> for method in METHODS:
...     np.vectorize(J_B)(y_squared, method)
array([-1.033016, -2.37644 , -3.229206, -3.458686, -2.874375, -1.659047,
       -1.138573, -0.837874, -0.639674, -0.50042 ])
array([-1.033015, -2.376438, -3.229206, -3.458678, -2.874366, -1.659047,
       -1.138573, -0.837874, -0.639674, -0.50042 ])
array([-21.416452, -15.652599, -10.752462,  -6.674717,  -3.414974,
        -1.766015,  -2.101377,  -3.512787,  -5.883763,  -9.171826])
array([8.274402, 6.852961, 5.324542, 3.629909, 1.592409, 1.592409,
       3.629909, 5.324542, 6.852961, 8.274402])
array([-4.879475, -5.299257, -4.535294, -2.666966, -0.360079, -0.472706,
       -0.498092, -0.4295  , -0.358932, -0.298332])
array([-4.242277, -4.526832, -4.016455, -2.663948, -0.616588, -0.506298,
       -0.513178, -0.43695 , -0.362923, -0.3006  ])

>>> for method in METHODS:
...     np.vectorize(J_B)(y_squared, method)
array([-1.033016, -2.37644 , -3.229206, -3.458686, -2.874375, -1.659047,
       -1.138573, -0.837874, -0.639674, -0.50042 ])
array([-1.033015, -2.376438, -3.229206, -3.458678, -2.874366, -1.659047,
       -1.138573, -0.837874, -0.639674, -0.50042 ])
array([-21.416452, -15.652599, -10.752462,  -6.674717,  -3.414974,
        -1.766015,  -2.101377,  -3.512787,  -5.883763,  -9.171826])
array([8.274402, 6.852961, 5.324542, 3.629909, 1.592409, 1.592409,
       3.629909, 5.324542, 6.852961, 8.274402])
array([-4.879475, -5.299257, -4.535294, -2.666966, -0.360079, -0.472706,
       -0.498092, -0.4295  , -0.358932, -0.298332])
array([-4.242277, -4.526832, -4.016455, -2.663948, -0.616588, -0.472706,
       -0.498092, -0.4295  , -0.358932, -0.298332])

Test their derivatives.

>>> for method in METHODS_DERIV:
...     np.vectorize(J_B)(y_squared, method=method, derivative=1)
...     np.vectorize(J_B)(y_squared, method=method, derivative=2)
array([-0.698182, -0.502894, -0.254838,  0.061861,  0.490156,  0.329235,
        0.170331,  0.107608,  0.073783,  0.053084])
array([ 0.077352,  0.099149,  0.126199,  0.163222,  0.23795 , -0.132829,
       -0.040093, -0.019849, -0.011616, -0.007431])
array([-0.698167, -0.502661, -0.254797,  0.062199,  0.490426,  0.329235,
        0.170331,  0.107608,  0.073783,  0.053084])
array([ 0.077264,  0.100877,  0.129686,  0.158361,  0.232503, -0.132827,
       -0.040094, -0.019848, -0.011616, -0.007431])

>>> for method in METHODS_DERIV:
...     np.vectorize(J_F)(y_squared, method=method, derivative=1)
...     np.vectorize(J_F)(y_squared, method=method, derivative=2)
array([ 2.164277,  0.287368, -0.2634  , -0.515941, -0.54348 , -0.262475,
       -0.153157, -0.101028, -0.070802, -0.051591])
array([ 0.261334, -0.350026, -0.166337, -0.062141,  0.047204,  0.078064,
        0.031681,  0.017306,  0.01063 ,  0.006992])
array([ 2.168401,  0.2869  , -0.263432, -0.515614, -0.543209, -0.262475,
       -0.153157, -0.101028, -0.070802, -0.051591])
array([ 0.766552, -0.396444, -0.162035, -0.064131,  0.040532,  0.078063,
        0.031681,  0.017306,  0.010631,  0.006992])

Test derivative against numerical derivatives. Look at positive y_squared, as
we know the numerical suffers from noise especially for negative arguments.

>>> y_squared = np.linspace(1, 10, 10)

>>> bessel = np.vectorize(J_F)(y_squared, derivative=1, method="bessel")
>>> approx = np.vectorize(J_F)(y_squared, derivative=1, method="approx")
>>> np.all(np.isclose(bessel, approx, rtol=1e-2))
True

>>> bessel = np.vectorize(J_F)(y_squared, derivative=2, method="bessel")
>>> approx = np.vectorize(J_F)(y_squared, derivative=2, method="approx")
>>> np.all(np.isclose(bessel, approx, rtol=1e-1))
True

>>> bessel = np.vectorize(J_B)(y_squared, derivative=1, method="bessel")
>>> approx = np.vectorize(J_B)(y_squared, derivative=1, method="approx")
>>> np.all(np.isclose(bessel, approx, rtol=1e-2))
True

>>> bessel = np.vectorize(J_B)(y_squared, derivative=2, method="bessel")
>>> approx = np.vectorize(J_B)(y_squared, derivative=2, method="approx")
>>> np.all(np.isclose(bessel, approx, rtol=1e-1))
True
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
