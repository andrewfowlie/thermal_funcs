"""
Plot thermal functions
======================
"""

import matplotlib.pyplot as plt
import numpy as np
import itertools
import time
from matplotlib import rc

from thermal_funcs import J_B, J_F

# LaTeX

rc('font',**{'family': 'serif', 'serif': ['Palatino']})
rc('text', usetex=True)

# Make plots

LINES = ["-", "--", "-.", ":", "-.", ":", "-.", ":"]
COLORS = ['Crimson', 'DarkSeaGreen', 'DodgerBlue', 'DarkSlateGray', 'Sienna', 'Gold', 'Violet']

def latex_float(float_):
    str_float = '{:.1E}'.format(float_)
    a, e = str_float.split('E')
    return r'${} \times 10^{{{}}}$'.format(a, int(e))

def plot_func(bosonic, methods, y_squared, name):

    fig = plt.figure(figsize=(8, 8))
    l_cycle = itertools.cycle(LINES)
    c_cycle = itertools.cycle(COLORS)

    for method in methods:

        func = J_B if bosonic else J_F
        t0 = time.time()
        J = np.array([func(y, method=method) for y in y_squared])
        t1 = time.time()

        av_time = (t1 - t0) / len(y_squared)

        label = r"{} --- mean time = {} s".format(method, latex_float(av_time))
        line = next(l_cycle)
        color = next(c_cycle)

        plt.plot(y_squared, J, label=label, ls=line, color=color, alpha=0.85, lw=3)


    if bosonic:
        plt.ylabel('$J_B(y^2)$')
    else:
        plt.ylabel('$J_F(y^2)$')

    plt.grid()
    plt.xlabel('$y^2$')
    leg = plt.legend(loc='best')
    leg.get_frame().set_alpha(0.75)
    plt.grid()
    plt.title('Thermal functions')
    plt.savefig(name)


plot_func(True, ['quad', 'bessel', 'approx', 'lim', 'zeta'], np.linspace(-1000., 0., 1000), 'J_B_neg.pdf')
plot_func(False, ['quad', 'bessel', 'approx', 'lim', 'zeta'], np.linspace(-1000., 0., 1000), 'J_F_neg.pdf')

plot_func(True, ['quad', 'bessel', 'approx', 'lim', 'zeta'], np.linspace(-100., 0., 1000), 'J_B_small_neg.pdf')
plot_func(False, ['quad', 'bessel', 'approx', 'lim', 'zeta'], np.linspace(-100., 0., 1000), 'J_F_small_neg.pdf')

plot_func(True, ['quad', 'bessel', 'approx', 'taylor', 'zeta'], np.linspace(-3., 3., 1000), 'J_B_small.pdf')
plot_func(False, ['quad', 'bessel', 'approx', 'taylor', 'zeta'], np.linspace(-3., 3., 1000), 'J_F_small.pdf')

plot_func(True, ['quad', 'bessel', 'approx'], np.linspace(0., 50., 1000), 'J_B_pos.pdf')
plot_func(False, ['quad', 'bessel', 'approx'], np.linspace(0., 50., 1000), 'J_F_pos.pdf')
