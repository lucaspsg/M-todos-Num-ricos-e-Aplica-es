import math
import numpy as np

# x = (e^-t).cos(t) -> x' = -t.x -(e^-t).y
# y = sin(t) -> y' = cos(t) = x/(e^-t)

class ImplicitEulerXY():

    def f(self, t, y):
        return np.array([-t*y[0] - math.exp(-t)*y[1], y[0]/math.exp(-t)])

    def SAM(self, t, dt, y):
        diff = 10
        i = 0
        root = y + dt * self.f(t, y)
        while i < 20 and diff < 0.0001:
            prev_root = root
            root = y + dt * self.f(t + dt, root)
            diff = np.linalg.norm(root - prev_root)
            i += 1
        return root

    def calculate_points(self, y0, t0, tf, ini_n, end_n, r):
        n = ini_n

        while n <= end_n:
            y = np.zeros((n, 2))
            y[0] = y0
            dt = (tf - t0)/(n - 1)
            t = t0 + dt
            index = 1

            while t <= tf:
                y[index] = self.SAM(t, dt, y[index - 1])
                t += dt
                index += 1


            print(y)
            n *= r

a = ImplicitEulerXY()

a.calculate_points([1, 0], 0, 5, 16, 16, 2)
