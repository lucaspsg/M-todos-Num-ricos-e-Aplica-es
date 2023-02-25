import matplotlib.pyplot as plt
import numpy as np

class ShowRkkStability():

    def newton(self, f, df, x0, tol=1e-5, max_iter=20):
        for i in range(max_iter):
            x = x0 - f(x0)/df(x0)
            if np.linalg.norm(x - x0) < tol:
                return x
            x0 = x

    def rk11_stable_region(self):
        def f(x):
            return x + 1

        def df(x):
            return 1

        return self.get_stable_region_points(f, df, [-2j, 1, 0j, -1])

    def rk22_stable_region(self):
        def f(x):
            return x**2/2 + x + 1

        def df(x):
            return x + 1

        return self.get_stable_region_points(f, df, [-2j, 2j, -1j - 1.5])

    def rk33_stable_region(self):
        def f(x):
            return x**3/6 + x**2/2 + x + 1

        def df(x):
            return x**2/2 + x + 1

        return self.get_stable_region_points(f, df, [-2, -1 +2j, -1 -2j, 2])

    def rk44_stable_region(self):
        def f(x):
            return x**4/24 + x**3/6 + x**2/2 + x + 1

        def df(x):
            return x**3/6 + x**2/2 + x + 1

        return self.get_stable_region_points(f, df, [-2 + 1j, -2 -1j, -1/2 + 5j/2, -1/2 -5j/2, -2j -2, -3 + 1j])

    def get_stable_region_points(self, f, df, x0):
        def F(x):
            thetas = np.linspace(-np.pi, np.pi, 100)
            return f(x) - np.exp(thetas * 1j)

        stable_region_points = np.array([])

        for ini_x in x0:
            stable_region_points = np.append(stable_region_points,
                                             self.newton(F, df, ini_x)
                                             )

        return stable_region_points

    def plot_stability_regions(self):
        plt.scatter(np.real(self.rk11_stable_region()),
                 np.imag(self.rk11_stable_region()),
                 label='RK11',
                 s=1
                 )

        plt.scatter(np.real(self.rk22_stable_region()),
                 np.imag(self.rk22_stable_region()),
                 label='RK22',
                 s=1
                 )

        plt.scatter(np.real(self.rk33_stable_region()),
                 np.imag(self.rk33_stable_region()),
                 label='RK33',
                 s=1
                 )

        plt.scatter(np.real(self.rk44_stable_region()),
                 np.imag(self.rk44_stable_region()),
                 label='RK44',
                 s=1
                 )

        plt.axis("scaled")
        plt.legend()
        plt.show()

a = ShowRkkStability()
a.plot_stability_regions()
