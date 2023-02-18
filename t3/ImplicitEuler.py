import math
import matplotlib.pyplot as plt
import numpy as np

# x = (e^-t).cos(t) -> x' = -e^(-t)(y + cos(t))
# y = sin(t) -> y' = cos(t) = x/(e^-t)

class ImplicitEulerXY():

    def f(self, t, y):
        return np.array([-math.exp(-t)*(y[1] + math.cos(t)),
                         y[0]/math.exp(-t)])

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
            t = np.zeros((n,))
            y[0] = y0
            t[0] = t0
            dt = (tf - t0)/n
            index = 1

            while index < n:
                t[index] = t[index - 1] + dt
                y[index] = self.SAM(t[index], dt, y[index - 1])
                index += 1


            if n <= 64:
                linestyle = ''
                if n == 16:
                    linestyle = 'solid'
                elif n == 32:
                    linestyle = 'dashed'
                else:
                    linestyle = 'dotted'

                plt.figure(0)
                plt.xlabel('t')
                plt.ylabel("y'(t)")
                plt.title("Aproximações de x(t) = (e^-t).cos(t)")

                plt.plot(t, y[:, 0], linestyle=linestyle, color='k',
                         label = 'n = ' + str(n))

                plt.figure(1)
                plt.xlabel('t')
                plt.ylabel("y(t)")
                plt.title('Aproximações de y(t) = sin(t)')

                plt.plot(t, y[:, 1], linestyle=linestyle, color='k',
                         label = 'n = ' + str(n))
                n *= r
                
        plt.figure(0)
        plt.legend(loc='best')
        plt.figure(1)
        plt.legend(loc='best')
        plt.show()

a = ImplicitEulerXY()

a.calculate_points([1, 0], 0, 3, 16, 64, 2)
