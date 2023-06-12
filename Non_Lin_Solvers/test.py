import numpy as np
def func(x):
    g = 9.806
    k = 0.00341
    f = np.log(np.cosh(x*np.sqrt(g*k)))- 50
    f2 = x
    return f

print(func(277.2210016846657))