"""
Time CosmoTransitions versus C++ implementation
===============================================
"""

from thermal_funcs import J_B, J_F
from cosmoTransitions.finiteT import _Jf_exact2, _Jb_exact2, Jf_high, Jb_high, Jf_low, Jb_low

import numpy as np
import timeit

n = 1000

# Timing for quadrature

y_squared = 10.
y = y_squared**0.5
setup = "from __main__ import y, y_squared, _Jb_exact2, _Jf_exact2, Jf_high, Jb_high, Jf_low, Jb_low, J_B, J_F"

python = min(timeit.repeat("_Jf_exact2(y_squared); _Jb_exact2(y_squared)", setup=setup, repeat=10, number=n)) / n
cpp = min(timeit.repeat("J_F(y_squared, method='quad'); J_B(y_squared, method='quad')", setup=setup, repeat=10, number=n)) / n

delta = [J_F(y_squared) + _Jf_exact2(y_squared), J_B(y_squared) - _Jb_exact2(y_squared)]

print "delta = {}".format(delta)
print "C++ min time  = {}. Python min time = {}".format(cpp, python)
# delta = [-1.3653467245688944e-11, 1.4502621326073495e-11]
# C++ min time  = 4.11880016327e-05. Python min time = 0.000731825113297

# Timing for Bessel

y_squared = 100.
y = y_squared**0.5
setup = "from __main__ import y, y_squared, _Jb_exact2, _Jf_exact2, Jf_high, Jb_high, Jf_low, Jb_low, J_B, J_F"

python = min(timeit.repeat("Jf_high(y); Jb_high(y)", setup=setup, repeat=10, number=n)) / n
cpp = min(timeit.repeat("J_F(y_squared, method='bessel', max_n=8); J_B(y_squared, method='bessel', max_n=8)", setup=setup, repeat=10, number=n)) / n

delta = [J_F(y_squared) + Jf_high(y), J_B(y_squared) - Jb_high(y)]

print "delta = {}".format(delta)
print "C++ min time  = {}. Python min time = {}".format(cpp, python)
# delta = [6.0715321659188248e-18, 5.2041704279304213e-18]
# C++ min time  = 5.73301315308e-06. Python min time = 2.7764081955e-05

# Timing for Taylor

y_squared = 0.1
y = y_squared**0.5
setup = "from __main__ import y, y_squared, _Jb_exact2, _Jf_exact2, Jf_high, Jb_high, Jf_low, Jb_low, J_B, J_F"

python = min(timeit.repeat("Jf_low(y); Jb_low(y)", setup=setup, repeat=10, number=n)) / n
cpp = min(timeit.repeat("J_F(y_squared, method='taylor', max_n=20); J_B(y_squared, method='taylor', max_n=20)", setup=setup, repeat=10, number=n)) / n

delta = [J_F(y_squared) + Jf_low(y), J_B(y_squared) - Jb_low(y)]

print "delta = {}".format(delta)
print "C++ min time  = {}. Python min time = {}".format(cpp, python)
# delta = [5.6349129806676501e-09, 1.8799717782513881e-08]
# C++ min time  = 5.5570602417e-06. Python min time = 3.29279899597e-05
