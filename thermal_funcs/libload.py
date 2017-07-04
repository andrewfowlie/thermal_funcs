"""
Load functions from a dynmaic library via the header file
=========================================================
"""

import ctypes
import numpy as np
import json
import sys
from collections import OrderedDict as odict
from future.utils import viewitems

import CppHeaderParser as cparse

if sys.version_info[0] < 3:
    import functools32 as functools
else:
    import functools

# Dictionary for string to ctype

TYPE_DICT = dict(double=ctypes.c_double,
                 int=ctypes.c_int,
                 bool=ctypes.c_bool)

# Decorators for scalar, vectorized and memoized functions

def asscalar(func):
    """
    Make function result a scalar, if possible, not a length 0 or 1 array.
    """
    @functools.wraps(func)
    def rfunc(*args, **kwargs):
        """
        Wrapped function that converts result from array to scalar
        if result is size == 1
        """
        res = func(*args, **kwargs)
        if res.size == 1:
            res = np.asscalar(res)
        return res

    return rfunc

decorate = lambda func: asscalar(np.vectorize(functools.lru_cache()(func)))

def make_arg_wrapper(name, default):
    """
    :returns: Part of a functions arguments with name = default
    """
    if default:
        strip = ''.join(default.split())
        convert = str(json.loads(strip))
        arg = '{}={}'.format(name, convert)
    else:
        arg = name
    return arg

def make_wrapper(default):
    """
    :returns: Wrapper for a function that adds default arguments
    """
    wrapped_arg_list = [make_arg_wrapper(*pair) for pair in viewitems(default)]
    wrapped_arg_str = ', '.join(wrapped_arg_list)
    arg_str = ', '.join(default.keys())
    wrapped_func = 'f_wrapped_gen = lambda f_load: lambda {}: f_load({})'.format(wrapped_arg_str, arg_str)
    return wrapped_func

def libload(lib_name, header_name):
    """
    Load functions from a header file/library.
    """
    header = cparse.CppHeader(header_name)
    lib = ctypes.CDLL(lib_name)

    f_dict = dict()

    for func in header.functions:

        f_name = func['name']
        f_type = func['rtnType']
        f_ctype = TYPE_DICT[f_type]

        p_types = [param['type'] for param in func['parameters']]
        p_ctypes = [TYPE_DICT[p_type] for p_type in p_types]
        p_default = odict((param['name'], param.get('default')) for param in func['parameters'])

        f_load = lib.__getattr__(f_name)
        f_load.argtypes = p_ctypes
        f_load.restype = f_ctype

        f_wrapped_str = make_wrapper(p_default)
        exec(f_wrapped_str, globals())  # NB that this exec call defines f_wrapped_gen
        f_wrapped = f_wrapped_gen(f_load)

        f_decorated = decorate(f_wrapped)
        f_decorated.__name__ = f_name
        f_decorated.__doc__ = "Python interface to {} loaded from {}".format(f_name, lib_name)

        f_dict[f_name] = f_decorated

    return f_dict
