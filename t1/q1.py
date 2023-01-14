import matplotlib.pyplot as plt
import numpy as np
import math

n = 100

t0 = -math.pi
tf = math.pi
t = 0

dt = (tf - t0)/n

linestyle = ''
label = ''

for m in range(1, 4):
    x = np.array([])
    y = np.array([])

    if m == 1:
        linestyle = 'solid'
        label = 'm = 1'
    elif m == 2:
        linestyle = 'dashed'
        label = 'm = 2'
    else:
        linestyle = 'dotted'
        label = 'm = 3'

    while dt*t + t0 <= tf:
        x = np.append(x, dt*t + t0)
        y = np.append(y, m*(dt*t + t0))
        t += 1

    y = np.cos(y)

    plt.plot(x, y, linestyle=linestyle, label=label)
    t = 0

plt.xlabel('t')
plt.ylabel('y')
plt.title('cos(mt)')
plt.legend(loc='best')

plt.show()
