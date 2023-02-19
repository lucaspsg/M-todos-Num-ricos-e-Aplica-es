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

class ImplicitTrapeziumXY():

    # função de discretização local que retorna um numpy array
    # onde o primeiro indice representa x' e o segundo, y'
    def f(self, t, y):
        a = 0.1
        b = 0.2
        r = 0.05
        m = 0.1
        return np.array([-y[0]*m + y[0]*y[1]*b,
                         y[1]*r - a*y[1]*y[0]])

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
        final_y = np.zeros((int(math.log(end_n/ini_n, r)) + 1, 2))

        iter_n = ini_n # inicialização do primeiro n de iteração
        iter_dt = (tf - t0)/iter_n #  inicialização do primeiro dt de iteração


        while iter_n <= end_n:
            y = np.zeros((iter_n + 1, 2)) # vetor com os valores de y em cada ponto
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
            if iter_n <= 16:
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

            iter_n *= r # atualizando n para a próxima iteração
            iter_dt /= r # atualizando dt para a próxima iteração


        # calculando os ps ao final de tudo
        for i in range(0, len(n) - 1):
            if i > 1:
                p[i] = math.log(float(error[i])/float(error[i+1]), r)
            else:
                p[i] = -1

        # printando tabela
        for i in range(0, len(n) - 1):
            print(str(n[i]) + " & " + self.formatNumber(dt[i]) + " & " +
                  self.formatNumber(error[i]) + " & " +  self.formatNumber(p[i]) +
                  " \\\\" + "\n")

        plt.figure(0)
        plt.legend(loc='best')
        plt.figure(1)
        plt.legend(loc='best')
        plt.show()

a = ImplicitTrapeziumXY()

a.calculate_points([3, 1], 0, 5, 4, 16384 * 2, 2)
