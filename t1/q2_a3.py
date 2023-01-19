# a equação a ser resolvida é a de Lotka-Volterra
# Pk' = -Pk.m + Pk.Hk.b
# H' = Hk.r - a.Hk.Pk
#
# a = 0.0005
# b = 0.1
# r = 0.04
# m = 0.2
# temos duas variaveis de estado -> Pk' e Hk''
# isso forma duas equações na forma de Cauchy
#
# resolveremos usando o Método de Euler
import numpy as np
import math
from decimal import Decimal

# função de discretização de Pk
def discretization_function_P(Pk, m, Hk, b):
    return -Pk*m + Pk*Hk*b;

# função de discretização de Hk
def discretization_function_H(Hk, r, a, Pk):
    return Hk*r - a*Hk*Pk

# função que calcula Pk+1
def calculate_Pkp1(Pk, m, Hk, b, Hn):
    return Pk + Hn * discretization_function_P(Pk, m, Hk, b)

# função que calcula Hk+1
def calculate_Hkp1(Hk, r, a, Pk, Hn):
    return Hk + Hn * discretization_function_H(Hk, r, a, Pk)

function_range_0 = 0 # começo do intervalo a ser avaliado
function_range_f = 3 # fim do intervalo a ser avaliado

n = 128 # número inicial de partições do intevalo
iter_n = n # variável auxiliar para iterar
f = 2 # fator multiplicativo de n para cada iteração

H0 = 3 # y(0)
P0 = 1 # y'(0)

a = 0.0005
b = 0.1
r = 0.04  
m = 0.2

Tk = 0 # t = 0

Yks = np.array([H0, P0]) # array com Pk e Hk

error = 0 # erro da iteração
errors = np.array([]) # array com os erros de cada iteração

p = '-' # o valor inicial de p é invalido pois são necessários 2 erros para 
        # computá-lo

index = 0 # variável auxiliar para acessar o erro na matrizes de erros

print("Tabela de Convergência Numérica")

while iter_n <= 128: # loop para testar diferentes quantidades de partições

    Hn = (function_range_f - function_range_0) / iter_n # dt

    while Tk < function_range_f: # enquanto Tk estiver dentro do intervalo
        # calcula Pk+1 e Hk+2 e os atualiza no array de variáveis de estado
        Hk = calculate_Hkp1(Yks[0], r, a, Yks[1], Hn) 
        Pk = calculate_Pkp1(Yks[1], m, Yks[0], b, Hn)
        Yks = [Hk, Pk]
        # atualiza Tk
        Tk = Tk + Hn

        print("Hk:", Hk, "Pk:", Pk)

    # reseta as variáveis para o próximo loop
    iter_n *= f
    Tk = 0
    Yks = [H0, P0]
    # atualiza o index
    index += 1




