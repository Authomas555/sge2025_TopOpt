import sys, os
sys.path.append(os.path.relpath("utils"))  # assume the root is src
from Density import interpolationLibrary, get_default_args

## test deritvative

import numpy as np

# check the interpolation library

def test_structure_interpolation_dictionnary():
    """ Check the interpolation library structure"""
    for keys in interpolationLibrary.keys():
        assert "function" in interpolationLibrary[keys]
        assert "derivative" in interpolationLibrary[keys]
        assert "condition" in interpolationLibrary[keys]
        assert "latex" in interpolationLibrary[keys]
    


def test_derivative_interpolation_dictionnary():
    """ Check the derivative in the interpolation library"""
    h = 1e-6
    for keys in interpolationLibrary.keys():
        x = np.random.rand(100)
        f = interpolationLibrary[keys]["function"]
        dicoval = get_default_args(f)
        for defaultArgs in dicoval.keys():
            dicoval[defaultArgs] = 10*np.random.rand(100)+1
        f = lambda x : interpolationLibrary[keys]["function"](x, **dicoval)
        df = lambda x : interpolationLibrary[keys]["derivative"](x, **dicoval)
        # second order finite difference
        assert np.allclose(df(x),  (f(x+h) - f(x-h)) / (2*h), atol = h, rtol = h), 'Derivative check failed for "' + keys + '"'



if __name__ == "__main__" : 
    print('testing utils/Density.py ...')
    test_structure_interpolation_dictionnary()
    test_derivative_interpolation_dictionnary()
    print('done !')