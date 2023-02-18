from decimal import Decimal
import math
import matplotlib.pyplot as plt
import numpy as np

# a equação a ser resolvida é a de Lotka-Volterra
# x' = -x.m + x.y.b
# y' = y.r - a.y.x
#
# H(0) = 3; P(0) = 1
#
# a solução será avaliada  o intervalo de 0 a 5

class ImplicitEulerXY():

    def f(self, t, y):
        a = 0.1
        b = 0.2
        r = 0.05
        m = 0.1
        return np.array([-y[0]*m + y[0]*y[1]*b,
                         y[1]*r - a*y[1]*y[0]])

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

    def formatNumber(self, n):
        if n == "-":
            return "----------"
        return str("{:.5E}".format(Decimal(n)))

    def calculate_points(self, y0, t0, tf, ini_n, end_n, r):
        iteration = 0
        n = ini_n
        error = np.zeros((int(math.log(end_n/ini_n, r)) + 1))
        p = "-"

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

            error[iteration] = 5

            if iteration > 0:
                p = math.log(error[iteration - 1]/error[iteration], r)

            print(str(n) + " & " + self.formatNumber(dt) + " & " +
                  self.formatNumber(error[iteration]) + " & " +
                  self.formatNumber(p) + " \\\\" + "\n")

            iteration += 1

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

a.calculate_points([3, 1], 0, 5, 16, 16384, 2)
