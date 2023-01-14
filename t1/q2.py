# o problema de solução inicial vai ser y'(t) = e^t.cos(t) com y(t0) = 0
# usaremos o método de Euler para resolvê-lo
# Yk+1 = Yk + dt.f(t, y(t)) = Yk + dt.y' = Yk + dt.(e^t.cos(t))
import numpy as np
import math

def exp(n, t):
    aux = n
    while(t > 1):
        n *= aux
        t -= 1

    return n

def discretization_function(t, y):
    pass
