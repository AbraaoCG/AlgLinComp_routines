yo = 1
to = 0
T = 2


def f(t,y):
    return -5*t*(y**2)/2

def rk2Ordem(yo,to,T):
    dt = 0.001
    N = T / dt
    yd1 = 0
    yi = 1
