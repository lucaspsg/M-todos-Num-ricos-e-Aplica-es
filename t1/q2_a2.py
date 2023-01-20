# a equação a ser resolvida é y'' + 4y' - 5y = 14 - 30t; y(0) = 3; y'(0) = 1
# para esse problema a solução é y(t) = 6t + e^(-5t) + 2
# usaremos o método de Euler para resolvê-lo
#
# temos duas variaveis de estado -> y' e y''
# y' = y1; y'' = y2
# isso forma duas equações na forma de Cauchy
#
# y1' = y2
# y2' = -4y1' + 5y1 + 14 - 30t = -4y2 + 5y1 + 14 - 30t
#
# Assim:
# Y2k+1 = Y2k + dt.r(Tk, Y2k, Y1k) = Y2k + dt.(-4Y2k + 5Y1k + 14 - 30Tk)
# Y1k+1 = Y1k + dt.r(Tk, Y2k, Y1k) = Y1k + dt.Y2k
import numpy as np
import math
from decimal import Decimal

# função de discretização de Y1
def discretization_function_y1(Y2k):
    return Y2k;

# função de discretização de Y2
def discretization_function_y2(Tk, Y1k, Y2k):
    return -4 * Y2k + 5 * Y1k + 14 - 30 * Tk

# função que calcula Y1k+1
def calculate_Y1kp1(Y1k, Y2k, Hn):
    return Y1k + Hn * discretization_function_y1(Y2k)

# função que calcula Y2k+1
def calculate_Y2kp1(Y1k, Y2k, Tk, Hn):
    return Y2k + Hn * discretization_function_y2(Tk, Y1k, Y2k)

# função de y
def y(t):
    return 6*t + math.exp(-5*t) + 2

# função de y'
def dy(t):
    return 6 - 5 * math.exp(-5*t)

#formata os números o output da tabela
def formatNumber(n):
    if n == "-":
        return "----------"
    return str("{:.5E}".format(Decimal(n)))

function_range_0 = 0 # começo do intervalo a ser avaliado
function_range_f = 3 # fim do intervalo a ser avaliado

n = 128 # número inicial de partições do intevalo
iter_n = n # variável auxiliar para iterar
r = 2 # fator multiplicativo de n para cada iteração

y0 = 3 # y(0)
dy0 = 1 # y'(0)

Tk = 0 # t = 0

Yks = np.array([y0, dy0]) # array com Y1k e Y2k

error = 0 # erro da iteração
errors = np.array([]) # array com os erros de cada iteração

p = '-' # o valor inicial de p é invalido pois são necessários 2 erros para 
        # computá-lo

index = 0 # variável auxiliar para acessar o erro na matrizes de erros

print("Tabela de Convergência Numérica")

while iter_n <= 16384: # loop para testar diferentes quantidades de partições

    Hn = (function_range_f - function_range_0) / iter_n # dt

    while Tk < function_range_f: # enquanto Tk estiver dentro do intervalo
        # calcula Y1k+1 e Y2k+2 e os atualiza no array de variáveis de estado
        Yk1 = calculate_Y1kp1(Yks[0], Yks[1], Hn) 
        Yk2 = calculate_Y2kp1(Yks[0], Yks[1], Tk, Hn)
        Yks = [Yk1, Yk2]
        # atualiza Tk
        Tk = Tk + Hn

    # calcula o erro usando a norma do máximo e o adiciona no array
    error1 = abs(Yks[0] - y(Tk))
    error2 = abs(Yks[1] - dy(Tk))
    error = max(error1, error2)

    errors = np.append(errors, error)

    # calcula p
    if iter_n > n:
        p = math.log(errors[index - 1]/error, r)

    # printa o resultado
    print(str(iter_n) + " & " + formatNumber(Hn) + " & " + formatNumber(error) + " & " +  formatNumber(p) + " \\\\" + "\n")

    # reseta as variáveis para o próximo loop
    iter_n *= r
    Tk = 0
    Yks = [y0, dy0]
    # atualiza o index
    index += 1




