import math
import numpy as np

# x = (e^-t).cos(t) -> x' = -t.x -(e^-t).y
# y = sin(t) -> y' = cos(t) = x/(e^-t)

class ImplicitEulerXY():

    def f(t, y):
        return np.array([[-t*y[0] - math.exp(-t)*y[1]], [y[0]/math.exp(-t)]])

    def SAM(t, dt, y):
        diff = 10
        i = 0
        root = y + dt * f(t, y)
        while i < 20 and diff < 0.0001:
            prev_root = root
            root = dt * f(t + dt, y)
            diff = np.linalg.norm(root - prev_root)
            i += 1
        return root

    def calculate_points(y0, t0, tf, ini_n, end_n, r):
        n = ini_n

        while n <= end_n:
            y = np.zeros(n, 2)
            y[0] = y0
            t = t0
            dt = (tf - t0)/n
            index = 1

            while t <= tf:
                y[index] = SAM(t, dt, y[index])
                t += dt
                index += 1

            y[index]


a = ImplicitEulerXY()

a.calculate_points([1, 2], 0, 5, 16, 256, 2)

