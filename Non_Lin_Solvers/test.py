import numpy as np
def func(x):
    g = 9.806
    k = 0.00341
    f = np.log10(np.cosh(x*np.sqrt(g*k)))- 50
    f2 = x
    return f

print(func(-633.387416601181))
