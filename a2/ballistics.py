from decimal import Decimal
import math
import matplotlib.pyplot as plt
import numpy as np

# a equação a ser resolvida é da posição de um corpo em um lancamento oblíquo levando em conta a resistência do ar
# x'' = -(k.x'.sqrt(x'^2+y'^2))/m
# y'' = -g -(k.y'.sqrt(x'^2 + y'^2))/m
#
# y0 = x -> f = y2
# y1 = y -> f = y3
# y2 = x' -> f = -u.y2.sqrt(y2^2 + y3^2)
# y3 = y' -> f = -g -u.y3.sqrt(y2^2 + y3^2)


class ImplicitTrapezium():

    # função de discretização local que retorna um numpy array
    # onde o primeiro indice representa x' e o segundo, y'
    def f(self, t, y):
        g = 9.8 # gravidade
        k = 0.82 * 1.28 * 0.00006 / 2 # c . p . A / 2
        m = 0.008 # massa
        u = k/m # mi
        return np.array([y[2], y[3], -u*y[2]*math.sqrt((y[2]**2 + y[3]**2)), -g - u*y[3]*math.sqrt((y[2]**2 + y[3]**2))])

    # implementação do método de aproximações sucessivas
    def SAM(self, t, dt, y):
        diff = 10 # diferença inicial entre cada iteração
        i = 0 # variável de contagem de iterações

        # root -> raiz a ser achada; o valor abaixo é o chute inicial
        root = y + (dt/2) * (self.f(t, y) + self.f(t + dt, y))

        # a cada loop faz Xn+1 = phi(Xn)
        while i < 10 and diff > 0.0001:
            prev_root = root
            root = y + (dt/2) * (self.f(t, y) + self.f(t + dt, root))
            diff = np.linalg.norm(root - prev_root)
            i += 1
        return root

    # implementacao do metodo spline cubica
    def interpolate_cubic_spline(self, x, y, target):
        n = len(x)
        if n != len(y):
            raise ValueError("Os arrays de entrada devem ter o mesmo tamanho.")
        
        # Encontre os índices do intervalo que contém o ponto de interesse.
        i = np.searchsorted(x, target)
        if i == 0:
            i = 1
        elif i == n:
            i = n - 1
        
        # Use os pontos vizinhos para construir a spline cúbica.
        x0, x1, x2, x3 = x[i-1:i+3]
        y0, y1, y2, y3 = y[i-1:i+3]
        
        h1 = x1 - x0
        h2 = x2 - x1
        h3 = x3 - x2
        
        d1 = (y1 - y0) / h1
        d2 = (y2 - y1) / h2
        d3 = (y3 - y2) / h3
        
        a = d1
        b = ((d2 - d1) / h2) - ((d1 - d3) / h3)
        c = (((d1 - d2) / h2) * h2 + ((d1 - 2*d2 + d3) / (h3**2)) * (h2**2)) / h2
        d = (((d1 - d2) / (h2**2)) * h2 + ((2*d1 - 5*d2 + 4*d3 - d3) / (h3**3)) * (h2**3)) / h2
        
        t = target - x1
        return y1 + a*t + b*(t**2) + c*(t**3) + d*(t**4)

    # formata os números para printagem da tabela
    def formatNumber(self, n):
        if n == "-":
            return "----------"
        return str("{:.5E}".format(Decimal(n)))

    def calculate_points(self, y0, t0, tf, ini_n, end_n, r):
        # conta as iterações (pra cada iteração, um n)
        iteration = 0

        # guarda o n de cada iteração
        n = np.zeros((int(math.log(end_n/ini_n, r)) + 1))

        # guarda o dt de cada iteração
        dt = np.zeros((int(math.log(end_n/ini_n, r)) + 1))

        # guarda o erro de cada iteração
        error = np.zeros((int(math.log(end_n/ini_n, r)) + 1))

        # guarda o p de cada iteração
        p = np.zeros((int(math.log(end_n/ini_n, r)) + 1))

        # guarda o último y de cada iteração
        final_y = np.zeros((int(math.log(end_n/ini_n, r)) + 1, len(y0)))

        iter_n = ini_n # inicialização do primeiro n de iteração
        iter_dt = (tf - t0)/iter_n #  inicialização do primeiro dt de iteração


        while iter_n <= end_n:
            y = np.zeros((iter_n + 1, len(y0))) # vetor com os valores de y em cada ponto
            t = np.zeros((iter_n + 1,)) #vetor com os valores de t em cada ponto
            y[0] = y0 # colocando os valores iniciais de y
            t[0] = t0 # colocando os valores iniciais de t
            index = 1 # index = 1, pois no index = 0 temos y0 e t0

            while index <= iter_n:
                t[index] = t[index - 1] + iter_dt # Tk+1 = Tk + dt

                # usando o Método das Aproximações Sucessivas para calcular
                # Yk+1
                y[index] = self.SAM(t[index - 1], iter_dt, y[index - 1])

                index += 1 # ir pro próximo ponto

            n[iteration] = iter_n # adicionando o n da iteração ao vetor
            dt[iteration] = iter_dt # adicionando o dt da iteração ao vetor
            final_y[iteration] = y[-1] # adicionando o ultimo y da iteração ao vetor

            # só calcular o erro se tiver y_final passado
            if iteration > 0:
                error[iteration] = abs(np.linalg.norm(final_y[iteration - 1]) -
                                   np.linalg.norm(final_y[iteration]))/3
            # -1 -> valor inválido
            else:
                error[iteration] = -1

            # vai pro próximo n
            iteration += 1

            # definindo o estilo de linha pra cada n na plotagem do gráfico
            if iter_n <= 64:
                linestyle = ''
                if iter_n == 4:
                    linestyle = 'solid'
                elif iter_n == 8:
                    linestyle = 'dashed'
                else:
                    linestyle = 'dotted'

                # plotagem do gráfico
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

                plt.figure(2)
                plt.xlabel('x')
                plt.ylabel("y")
                plt.title('Aproximações de y(t)')

                plt.plot(y[:, 0], y[:, 1], linestyle=linestyle, color='r',
                         label = 'n = ' + str(iter_n))

            iter_n *= r # atualizando n para a próxima iteração
            iter_dt /= r # atualizando dt para a próxima iteração


        # calculando os ps ao final de tudo
        for i in range(0, len(n) - 1):
            if i > 1:
                p[i] = math.log(float(error[i])/float(error[i+1]), r)
            else:
                p[i] = -1

         #printando tabela
        for i in range(0, len(n) - 1):
            print(str(n[i]) + " & " + self.formatNumber(dt[i]) + " & " +
                self.formatNumber(error[i]) + " & " +  self.formatNumber(p[i]) +
                   " \\\\" + "\n")

        plt.figure(0)
        plt.legend(loc='best')
        plt.figure(1)
        plt.legend(loc='best')
        plt.show()
        print(self.interpolate_cubic_spline(y[:,0], y[:,1], 42.5))


a = ImplicitTrapezium()

a.calculate_points([0, 1, 50, 0], 0, 1, 4, 16384 * 2, 2)
