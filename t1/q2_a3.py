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

#formata os números o output da tabela
def formatNumber(n):
    if n == "-":
        return "----------"
    return str("{:.5E}".format(Decimal(n)))

function_range_0 = 0 # começo do intervalo a ser avaliado
function_range_f = 3 # fim do intervalo a ser avaliado

n = 128 # número inicial de partições do intevalo
Ns = np.array([]) # array com os Ns
iter_n = n # variável auxiliar para iterar
f = 2 # fator multiplicativo de n para cada iteração

Hns = np.array([]) # array com os Hns

H0 = 3 # y(0)
P0 = 1 # y'(0)

a = 0.0005
b = 0.1
r = 0.04  
m = 0.2

Tk = 0 # t = 0

Hs = np.array([]) # array com os Hks
Ps = np.array([]) # array com os Pks

error = 0 # variável que guarda o erro de norma euclidiana
errors = np.array([]) # array com os erros

p = '-'
ps = np.array(['-'])

index = 0 # variável auxiliar para acessar o erro na matrizes de erros

print("Tabela de Convergência Numérica")

while iter_n <= 16384: # loop para testar diferentes quantidades de partições
    Ns = np.append(Ns, iter_n)

    Hk = H0
    Pk = P0

    Hn = (function_range_f - function_range_0) / iter_n # dt
    Hns = np.append(Hns, Hn)

    while Tk < function_range_f: # enquanto Tk estiver dentro do intervalo
        # calcula Pk+1 e Hk+2 e os atualiza no array de variáveis de estado
        Hk = calculate_Hkp1(Hk, r, a, Pk, Hn)
        Pk = calculate_Pkp1(Pk, m, Hk, b, Hn)


        # atualiza Tk
        Tk = Tk + Hn

    Hs = np.append(Hs, Hk)
    Ps = np.append(Ps, Pk)

    if index > 0:
        # calcula o erro usando a norma euclidiana e a adiciona no array
        errorHkm1 = -(Hs[-1] - Hs[-2])
        errorPkm1 = -(Ps[-1] - Ps[-2])
        error = errorHkm1**2 + errorPkm1**2
        errors = np.append(errors, math.sqrt(error))

    # reseta as variáveis para o próximo loop
    iter_n *= f
    Tk = 0
    Hk = H0
    Pk = P0
    # atualiza o index
    index += 1


for i in range(0, len(Hs) - 2):
    p = math.log(errors[i]/errors[i+1], f)
    ps = np.append(ps, p)

errors = np.append(errors, '-')
ps = np.append(ps, '-')

for i in range(0, len(Hs)):
    print(str(Ns[i]) + " & " + formatNumber(Hns[i]) + " & " + formatNumber(errors[i]) + " & " +  formatNumber(ps[i]) + " \\\\" + "\n")
