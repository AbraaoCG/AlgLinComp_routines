import numpy as np
def func(x):
    g = 9.806
    k = 0.00341
    f = np.log10(np.cosh(x*np.sqrt(g*k)))- 50
    f2 = ( ( 4 * np.cos(x) ) - ( np.e**(2*x) ) )
    return f2

print(func(0))