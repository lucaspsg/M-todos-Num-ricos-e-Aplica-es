# x² - 3x + 2 = f(x)
# x1 = 1 ; x2 = 2 ;
# x² - 3x + 2 = 0 -> x(x - 3) + 2 = 0 -> x = 2/(x - 3) 
# I1 = [0.5 ; 1.5]
# I2 = [1.5 ; 2.5]

def ponto_fixo(x):
    x_prev = -999
    while abs(x - x_prev) > 0.00001:
        print(x)
        x_prev = x
        x = (3*x - 2)/x

ponto_fixo(2.5)

#fi1 = x² - 2x + 2 -> fi1' = 2x - 2
#fi2 = (3x - 2)/x -> fi2' = (3*x - 3x*ln(x))/x**2
#fi3 =  -2/(x - 3) -> -2*ln(x - 3) 
#fi4 = sqrt(3x - 2) -> fi4' = (1/2)*(3x - 2)**(-1/2)*3
#fi5
