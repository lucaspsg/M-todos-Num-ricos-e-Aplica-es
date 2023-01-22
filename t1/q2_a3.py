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

Hk = H0
Pk = P0

# coeficientes do Lotka-Volterra
a = 0.0005
b = 0.1
r = 0.04
m = 0.2

Tk = 0 # t = 0

Hs = np.array([]) # array com os Hks
Ps = np.array([]) # array com os Pks

error = 0 # variável que guarda o erro de norma euclidiana
errors = np.array(['-']) # array com os erros

p = '-' # variável que contém p
ps = np.array(['-', '-']) # array que contém os diferentes ps

index = 0 # variável auxiliar para acessar o erro na matrizes de erros

print("Tabela de Convergência Numérica")

while iter_n <= 16384: # loop para testar diferentes quantidades de partições
    Ns = np.append(Ns, iter_n) # adicionando n ao array

    Hn = (function_range_f - function_range_0) / iter_n # dt
    Hns = np.append(Hns, Hn) # adicionando dt ao array

    while Tk < function_range_f: # enquanto Tk estiver dentro do intervalo
        # calcula Pk+1 e Hk+2 e os atualiza no array de variáveis de estado
        Hk = calculate_Hkp1(Hk, r, a, Pk, Hn)
        Pk = calculate_Pkp1(Pk, m, Hk, b, Hn)


        # atualiza Tk
        Tk = Tk + Hn

    # salva os valores de yk no último k
    Hs = np.append(Hs, Hk)
    Ps = np.append(Ps, Pk)

    if index > 0:
        # calcula o erro usando a norma euclidiana e a adiciona no array
        errorHkm1 = -(Hs[-2] - Hs[-1])
        errorPkm1 = -(Ps[-2] - Ps[-1])
        error = math.sqrt(errorHkm1**2 + errorPkm1**2)
        errors = np.append(errors, error)

    # reseta as variáveis para o próximo loop
    iter_n *= f
    Tk = 0
    Hk = H0
    Pk = P0

    # atualiza o index
    index += 1


# calcula p agora que temos todos os erros
for i in range(1, len(Hs) - 1):
    p = math.log(float(errors[i])/float(errors[i+1]), f)
    ps = np.append(ps, p)

#printa a tabela de convergência
for i in range(0, len(Hs)):
    print(str(Ns[i]) + " & " + formatNumber(Hns[i]) + " & " + formatNumber(errors[i]) + " & " +  formatNumber(ps[i]) + " \\\\" + "\n")
