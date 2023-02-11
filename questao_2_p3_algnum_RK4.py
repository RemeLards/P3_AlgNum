from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
import math as m

# O resultado aproximado é obtido apenas por janelas pequenas de chutes  z <= |0.4| (mais que isso RK4 explode e Y se torna infinito muito rápido),
# Creio eu que isso é devido ao problema ser não linear, e portanto os passos "h" deveriam ser menores para resolver esse problema
# Como isso não é possivel para a prova, apenas fui chutando os valores (após obter uma noção através do algoritmo escrito)
# Encontrei um valor de "z" bem próximo, portanto usaremos valores bem próximos a ele para os dois chutes "zi" e "zi2"
# Valor de "z" encontrado no chute = -0.1763
# Valor de "z" encontrado pelo algoritmo chutando "zi" e "zi2" quando z <= |0.4| = -0.16490798660203138
# No caso acima, eu usei "zi" = 0.05 e "zi2" = "-zi"

#Equações desacopladas

# O valor da função z é dado pelo RK4 da função "dzdx_func", não temos uma descrição analítica para a função z
# Porém temos os valores da função por causa do RK4 da sua "derivada"
# Faço essa função assim pois o valor está armazenado no ultimo índice do vetor
# E também porque eu REAPROVEITO a estrutura do RK4, não precisando assim reescreve-lo

# Resumo, o RK4 de "dzdx_func" nos retorna o valor aproximado de "z"
# O RK4 de "z_func", nos retorna "y", uma vez que desacoplamos fazendo "dy/dx = z"
def z_func1(x=None, y=None):
    return dzdx_approx_vector[-1]

def z_func2(x=None, y=None):
    return dzdx_approx_vector2[-1]

def z_func3(x=None, y=None):
    return dzdx_approx_vector3[-1]

def z_func4(x=None, y=None):
    return dzdx_approx_vector4[-1]


def dzdx_func1(x ,y):
    return (y**3 - y*dzdx_approx_vector[-1])

def dzdx_func2(x,y):
    return (y**3 - y*dzdx_approx_vector2[-1])

def dzdx_func3(x,y):
    return (y**3 - y*dzdx_approx_vector3[-1])

def dzdx_func4(x,y):
    return (y**3 - y*dzdx_approx_vector4[-1])


def rungeKutta(x0, y0, h,func):
    
    y = y0 # Na primeira iteração y0 = y
   
    #Aplicando RK4

    k1 = h * func(x0, y)
    k2 = h * func(x0 + 0.5 * h, y + 0.5 * k1)
    k3 = h * func(x0 + 0.5 * h, y + 0.5 * k2)
    k4 = h * func(x0 + h, y + k3)

    # Update next value of y
    y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4)

    return y

# Fazendo o Primeiro chute
# shooting method
xi = 1 # Boundary Conditions
xf = 2 # Boundary Conditions 
h = 0.1
yi = 1/2 # Boundary Conditions
yf = 1/3 # Boundary Conditions 

zi = 0.4 # starting value guess 
x0 = xi # xi but iterative


#Vetores Para Armazenar todos os valores, e assim conseguir analizar o erro depois (e plotar gráfico)
dzdx_approx_vector = [zi]
z_approx_vector = [yi]

# Count number of iterations using step size or
# step height h
n = (int)((xf - xi)/h)


#Faço uma interação, para depois usar os próprios vetores para iterar (acho mais elegante)
dzdx_approx_vector.append(rungeKutta(x0, zi, h, dzdx_func1))

z_approx_vector.append(rungeKutta(x0, yi, h, z_func1))

x0 += h

i=1
while(i != n):
    
    dzdx_approx_vector.append(rungeKutta(x0, dzdx_approx_vector[i], h, dzdx_func1))
    
    z_approx_vector.append(rungeKutta(x0, z_approx_vector[i], h, z_func1))
    
    x0+= h # Update next value of x

    i+=1

#Vetor de x apenas para plotar
x_vector = [xi]
for i in range(n):
    x_vector.append(x_vector[i] + h)

print("\n")
print("Vetor X = " + str(x_vector))
print("Vetor Y_Chute1 = "  + str(z_approx_vector))

#Fazendo o Segundo chute

zi2 = -0.4
x0 = 1

dzdx_approx_vector2 = [zi2]
z_approx_vector2 = [yi]

#Faço uma interação, para depois usar os próprios vetores para iterar (acho mais elegante)
dzdx_approx_vector2.append(rungeKutta(x0, zi2, h, dzdx_func2))

z_approx_vector2.append(rungeKutta(x0, yi, h, z_func2))

x0 += h
i=1
while(i != n):
    
    dzdx_approx_vector2.append(rungeKutta(x0, dzdx_approx_vector2[i], h, dzdx_func2))
    
    z_approx_vector2.append(rungeKutta(x0, z_approx_vector2[i], h, z_func2))
    
    x0+= h # Update next value of x

    i+=1

print("\n")
print("Vetor X = " + str(x_vector))
print("Vetor Y_Chute2 = "  + str(z_approx_vector2))

# por interpolação temos que:

z_interpolado = zi + ( (zi2-zi) / (z_approx_vector2[-1] - z_approx_vector[-1]) ) * ( yf - z_approx_vector[-1] )
print("\n")
print("Z inicial real = " + str(z_interpolado))
print("\n")
x0 = 1

#Fazendo o gráfico aproximado 
dzdx_approx_vector3 = [z_interpolado]
z_approx_vector3 = [yi]

#Faço uma interação, para depois usar os próprios vetores para iterar (acho mais elegante)
dzdx_approx_vector3.append(rungeKutta(x0, z_interpolado, h, dzdx_func3))

z_approx_vector3.append(rungeKutta(x0, yi, h, z_func3))

x0 += h
i=1
while(i != n):
    
    dzdx_approx_vector3.append(rungeKutta(x0, dzdx_approx_vector3[i], h, dzdx_func3))
    
    z_approx_vector3.append(rungeKutta(x0, z_approx_vector3[i], h, z_func3))
    
    x0+= h # Update next value of x

    i+=1

print("\n")
print("Vetor X = " + str(x_vector))
print("Vetor Y_Interpolacao = "  + str(z_approx_vector3))

# Fazendo uso do chute "calculado"

zi3 = -0.1763
x0 = 1

dzdx_approx_vector4 = [zi3]
z_approx_vector4 = [yi]

#Faço uma interação, para depois usar os próprios vetores para iterar (acho mais elegante)
dzdx_approx_vector4.append(rungeKutta(x0, zi3, h, dzdx_func4))

z_approx_vector4.append(rungeKutta(x0, yi, h, z_func4))

x0 += h
i=1
while(i != n):
    
    dzdx_approx_vector4.append(rungeKutta(x0, dzdx_approx_vector4[i], h, dzdx_func4))
    
    z_approx_vector4.append(rungeKutta(x0, z_approx_vector4[i], h, z_func4))
    
    x0+= h # Update next value of x

    i+=1

print("\n")
print("Vetor X = " + str(x_vector))
print("Vetor Y_Chute3 = "  + str(z_approx_vector4))

y_real_sol_vec = []
for j in x_vector:
    y_real_sol_vec.append((j+1)**(-1))


error_vector = [0]
for k in range(n-1):
    error_vector.append( (1 - z_approx_vector4[k+1]/y_real_sol_vec[k+1]) * 100 )

print("\n")
print("Vetor Erro = "  + str(error_vector))
print("\n")
print("Erro Maximo = " + str(max(error_vector)) + "%")

erro_medio = 0
for l in range(n):
    erro_medio += error_vector[k]
erro_medio /= 10

print("Erro Medio = " + str(erro_medio) + "%")



figure, (axis1,axis2,axis3,axis4,axis5) = plt.subplots(1,5)

axis1.plot(x_vector,z_approx_vector)
axis1.set_title("Primeiro Chute")

axis2.plot(x_vector,z_approx_vector2)
axis2.set_title("Segundo Chute")

axis3.plot(x_vector,z_approx_vector3)
axis3.set_title("Resultado Interpolado")

axis4.plot(x_vector,z_approx_vector4)
axis4.set_title("Resultado Aproximado")

axis5.plot(x_vector,y_real_sol_vec)
axis5.set_title("Resultado Exato")

plt.show()
