import matplotlib.pyplot as plt
import numpy as np
import math

n = 100 # número de partições

# definindo os intervalos de partição
t0 = -math.pi
tf = math.pi
t = 0

dt = (tf - t0)/n # dt

# configuração dos gráficos
linestyle = ''
label = ''

# testando para os Ms pedidos
for m in range(1, 4):
    # arrays de valores para x e y
    x = np.array([])
    y = np.array([])

    # estilos de gráfico para cada m
    if m == 1:
        linestyle = 'solid'
        label = 'm = 1'
    elif m == 2:
        linestyle = 'dashed'
        label = 'm = 2'
    else:
        linestyle = 'dotted'
        label = 'm = 3'

    # enquanto t estiver no intervalo, fazer Yk = Yk+1
    while dt*t + t0 <= tf:
        x = np.append(x, dt*t + t0)
        y = np.append(y, m*(dt*t + t0))
        # atualizar t
        t += 1

    # plotando a função cos de fato
    y = np.cos(y)

    plt.plot(x, y, linestyle=linestyle, label=label, color='k')
    t = 0

plt.xlabel('t')
plt.ylabel('y')
plt.title('cos(mt)')
plt.legend(loc='best')

plt.show()
