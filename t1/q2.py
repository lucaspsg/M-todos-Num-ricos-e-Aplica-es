# o problema de solução inicial vai ser y'(t) = y - e^t.sin(t) com y(t0) = 0
# para esse problema, a solução é y(t) = e^t.cos(t)
# usaremos o método de Euler para resolvê-lo
# Yk+1 = Yk + Hn.f(t, y(t)) = Yk + Hn.y' = Yk + Hn.(Yk - e^Tk.sin(Tk))
import numpy as np
import math
import matplotlib.pyplot as plt
from decimal import Decimal

def discretization_function(Tk, Yk):
    return Yk - math.exp(Tk) * math.sin(Tk)

def calculate_Ykp1(Yk, Tk, Hn):
    return Yk + Hn * discretization_function(Tk, Yk)

def formatNumber(n):
    if n == "-":
        return "-"
    return str("{:.5E}".format(Decimal(n)))

function_range_0 = 0
function_range_f = 5

n = 128
ns = np.array([])
f = 2

Hns = np.array([])

y0 = 1

Tk = 0
Yk = y0

error = 0
errors = np.array([])

p = 0
ps = np.array(['-'])

index = 0

while n <= 16384:

    Hn = (function_range_f - function_range_0) / n
    t = np.array([])
    y = np.array([])

    while Tk < function_range_f:
        t = np.append(t, Tk)
        y = np.append(y, Yk)
        Yk = calculate_Ykp1(Yk, Tk, Hn)
        Tk = Tk + Hn

    error = abs(Yk - math.exp(Tk) * math.cos(Tk))

    ns = np.append(ns, n)
    Hns = np.append(Hns, Hn)
    errors = np.append(errors, error)

    if n > 128:
        ps = np.append(ps, math.log(errors[index - 1]/error, f))

    n *= f
    Tk = 0
    Yk = y0
    index += 1


f = open('./latex_output.o', "w")

for i in range(0, len(ns)):
    f.write(str(ns[i]) + " & " + formatNumber(Hns[i]) + " & " + formatNumber(errors[i]) + " & " +  formatNumber(ps[i]) + " \\\\" + "\n")

