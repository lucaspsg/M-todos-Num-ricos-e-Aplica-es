from decimal import Decimal
import math
import matplotlib.pyplot as plt
import numpy as np

# x = (e^-t).cos(t) -> x' = -e^(-t)(sin(t) - cos(t) + 2*cos(t)) -> x' = -2x - (e^-t) * y
# y = sin(t) - cos(t) -> y' = 2*cos(t) + sin(t) - cos(t) = 2*x/(e^-t) + y

class ImplicitTrapeziumXY():

    def f(self, t, y):
        return np.array([-2*y[0] - math.exp(-t)*y[1],
                         2*y[0]/math.exp(-t) + y[1]])

    def solution(self, t):
        return [math.exp(-t)*math.cos(t), math.sin(t) - math.cos(t)]

    def SAM(self, t, dt, y):
        diff = 10
        i = 0
        root = y + (dt/2) * (self.f(t, y) + self.f(t + dt, y))
        while i < 10 and diff > 0.0001:
            prev_root = root
            root = y + (dt/2) * (self.f(t, y) + self.f(t + dt, root))
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
            y = np.zeros((n + 1, 2))
            t = np.zeros((n + 1,))
            y[0] = y0
            t[0] = t0
            dt = (tf - t0)/n
            index = 1

            while index <= n:
                t[index] = t[index - 1] + dt
                y[index] = self.SAM(t[index - 1], dt, y[index - 1])
                #y[index] = y[index - 1] + dt * self.f(t[index - 1], y[index - 1])
                index += 1

            error[iteration] = np.linalg.norm(y[index - 1]
                                               - self.solution(t[index - 1]))

            if iteration > 0:
                p = math.log(error[iteration - 1]/error[iteration], r)

            print(str(n) + " & " + self.formatNumber(dt) + " & " +
                  self.formatNumber(error[iteration]) + " & " +
                  self.formatNumber(p) + " \\\\" + "\n")

            iteration += 1

            if n <= 16:
                linestyle = ''
                if n == 4:
                    linestyle = 'solid'
                elif n == 8:
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
                plt.title('Aproximações de y(t) = sin(t) - cos(t)')

                plt.plot(t, y[:, 1], linestyle=linestyle, color='k',
                         label = 'n = ' + str(n))
            n *= r

        plt.figure(0)
        plt.legend(loc='best')
        plt.figure(1)
        plt.legend(loc='best')
        plt.show()

a = ImplicitTrapeziumXY()

a.calculate_points([1, -1], 0, 3, 4, 16384, 2)
