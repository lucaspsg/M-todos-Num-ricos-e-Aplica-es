from decimal import Decimal
import math
import matplotlib.pyplot as plt
import numpy as np

# x = (e^-t).cos(t) -> x' = -e^(-t)(sin(t) - cos(t) + 2*cos(t)) -> x' = -2x - (e^-t) * y
# y = sin(t) - cos(t) -> y' = 2*cos(t) + sin(t) - cos(t) = 2*x/(e^-t) + y

class ImplicitTrapeziumXY():

    # função de discretização local que retorna um numpy array
    # onde o primeiro indice representa x' e o segundo, y' 
    def f(self, t, y):
        return np.array([-2*y[0] - math.exp(-t)*y[1],
                         2*y[0]/math.exp(-t) + y[1]])

    # funcao que retorna a solucao conhecida -> primeiro índice
    # representa a solução de x e o segundo, de y
    def solution(self, t):
        return [math.exp(-t)*math.cos(t), math.sin(t) - math.cos(t)]

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

    # execução do método do trapézio 
    def calculate_points(self, y0, t0, tf, ini_n, end_n, r):
        iteration = 0 # conta as iterações (pra cada iteração um n)
        n = ini_n # n inicial

        # erro de cada iteração
        error = np.zeros((int(math.log(end_n/ini_n, r)) + 1))

        p = "-" # p inicial

        while n <= end_n:
            y = np.zeros((n + 1, 2)) # vetor com os valores de y em cada ponto
            t = np.zeros((n + 1,)) # vetor com os valores de t em cada ponto
            y[0] = y0 # colocando os valores iniciais de y
            t[0] = t0 # colocando os valores iniciais de t
            dt = (tf - t0)/n # dt = h
            index = 1 # index = 1, pois no index = 0 temos y0 e t0

            while index <= n:
                t[index] = t[index - 1] + dt # Tk+1 = Tk + dt

                # usando o Método das Aproximações Sucessivas para calcular
                # Yk+1
                y[index] = self.SAM(t[index - 1], dt, y[index - 1])

                index += 1 # ir pro próximo ponto

            # erro = normal da diferença entre o resultado numérico e o 
            # resultado real
            error[iteration] = np.linalg.norm(y[index - 1]
                                               - self.solution(t[index - 1]))

            # só calcular p quando tiver erro atual e anterior
            if iteration > 0:
                p = math.log(error[iteration - 1]/error[iteration], r)

            # printa tabela
            print(str(n) + " & " + self.formatNumber(dt) + " & " +
                  self.formatNumber(error[iteration]) + " & " +
                  self.formatNumber(p) + " \\\\" + "\n")

            # vai pro próximo n
            iteration += 1

            # definindo o estilo de linha pra cada n na plotagem do gráfico
            if n <= 16:
                linestyle = ''
                if n == 4:
                    linestyle = 'solid'
                elif n == 8:
                    linestyle = 'dashed'
                else:
                    linestyle = 'dotted'

                # plotagem do gráfico
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

            n *= r # atualizando n para a próxima iteração

        plt.figure(0)
        plt.legend(loc='best')
        plt.figure(1)
        plt.legend(loc='best')
        plt.show()

a = ImplicitTrapeziumXY()

a.calculate_points([1, -1], 0, 3, 4, 16384, 2)
