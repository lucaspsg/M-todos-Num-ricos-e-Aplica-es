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
        if n == -1:
            return "----------"
        return str("{:.5E}".format(Decimal(n)))

    def calculate_points(self, y0, t0, tf, ini_n, end_n, r):
        iteration = 0
        n = np.zeros((int(math.log(end_n/ini_n, r)) + 1))
        dt = np.zeros((int(math.log(end_n/ini_n, r)) + 1))
        error = np.zeros((int(math.log(end_n/ini_n, r)) + 1))
        p = np.zeros((int(math.log(end_n/ini_n, r)) + 1))

        iter_n = ini_n
        iter_dt = (tf - t0)/iter_n

        while iter_n <= end_n:
            y = np.zeros((iter_n, 2))
            t = np.zeros((iter_n,))
            y[0] = y0
            t[0] = t0
            index = 1

            while index < iter_n:
                t[index] = t[index - 1] + iter_dt
                y[index] = self.SAM(t[index], iter_dt, y[index - 1])
                index += 1

            n[iteration] = iter_n
            dt[iteration] = iter_dt

            if iteration > 0:
                error[iteration] = np.linalg.norm(y[iteration - 1] -
                                                  y[iteration])
            else:
                error[iteration] = -1

            iteration += 1

            if iter_n <= 64:
                linestyle = ''
                if iter_n == 16:
                    linestyle = 'solid'
                elif iter_n == 32:
                    linestyle = 'dashed'
                else:
                    linestyle = 'dotted'

                plt.figure(0)
                plt.xlabel('t')
                plt.ylabel("y'(t)")
                plt.title("Aproximações de x(t)")

                plt.plot(t, y[:, 0], linestyle=linestyle, color='k',
                         label = 'n = ' + str(iter_n))

                plt.figure(1)
                plt.xlabel('t')
                plt.ylabel("y(t)")
                plt.title('Aproximações de y(t)')

                plt.plot(t, y[:, 1], linestyle=linestyle, color='k',
                         label = 'n = ' + str(iter_n))
            iter_n *= r
            iter_dt /= r

        for i in range(0, len(n) - 1):
            if i > 1:
                p[i] = math.log(float(error[i])/float(error[i+1]), r)
            else:
                p[i] = -1

        for i in range(0, len(n) - 1):
            print(str(n[i]) + " & " + self.formatNumber(dt[i]) + " & " +
                  self.formatNumber(error[i]) + " & " +  self.formatNumber(p[i]) +
                  " \\\\" + "\n")

        plt.figure(0)
        plt.legend(loc='best')
        plt.figure(1)
        plt.legend(loc='best')
        plt.show()

a = ImplicitEulerXY()

a.calculate_points([3, 1], 0, 5, 16, 16384 * 2, 2)
