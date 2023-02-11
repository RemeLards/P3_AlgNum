from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt
import math as m


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

def z_func_real(x=None, y=None):
    return dzdx_approx_vector3[-1]

def dzdx_func(x,y):
    return (4*(y-x))

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
xi = 0 # Boundary Conditions
xf = 1 # Boundary Conditions 
h = 0.1
yi = 0 # Boundary Conditions
yf = 2 # Boundary Conditions 

zi = -40  # starting value guess
x0 = 0 # x but iterative

#zi = 1 -> yf = 12.987136071093737
#zi2 = 0.2 -> yf = -0.01247573807291627
#calculei e deu zi = 0.32384835901969455, precisão não muito boa, aproximei bastante

#Vetores Para Armazenar todos os valores, e assim conseguir analizar o erro depois (e plotar gráfico)
dzdx_approx_vector = [zi]
z_approx_vector = [yi]

# Count number of iterations using step size or
# step height h
n = (int)((xf - x0)/h)


#Faço uma interação, para depois usar os próprios vetores para iterar (acho mais elegante)
dzdx_approx_vector.append(rungeKutta(x0, zi, h, dzdx_func))

z_approx_vector.append(rungeKutta(x0, yi, h, z_func1))

x0+= h

i=1
while(i != n):
    
    dzdx_approx_vector.append(rungeKutta(x0, dzdx_approx_vector[i], h, dzdx_func))
    
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

zi2 = 25
x0 = 0

dzdx_approx_vector2 = [zi2]
z_approx_vector2 = [yi]

#Faço uma interação, para depois usar os próprios vetores para iterar (acho mais elegante)
dzdx_approx_vector2.append(rungeKutta(x0, zi2, h, dzdx_func))

z_approx_vector2.append(rungeKutta(x0, yi, h, z_func2))

x0+= h
i=1
while(i != n):
    
    dzdx_approx_vector2.append(rungeKutta(x0, dzdx_approx_vector2[i], h, dzdx_func))
    
    z_approx_vector2.append(rungeKutta(x0, z_approx_vector2[i], h, z_func2))
    
    x0+= h # Update next value of x

    i+=1

print("\n")
print("Vetor X = " + str(x_vector))
print("Vetor Y_Chute2 = "  + str(z_approx_vector2))

# por interpolação temos que:

z_real = zi + ( (zi2-zi) / (z_approx_vector2[-1] - z_approx_vector[-1]) ) * ( yf - z_approx_vector[-1] )
print("\n")
print("Z inicial real = " + str(z_real))
print("\n")

x0 = 0

#Fazendo o gráfico aproximado 
dzdx_approx_vector3 = [z_real]
z_approx_vector3 = [yi]

#Faço uma interação, para depois usar os próprios vetores para iterar (acho mais elegante)
dzdx_approx_vector3.append(rungeKutta(x0, z_real, h, dzdx_func))

z_approx_vector3.append(rungeKutta(x0, yi, h, z_func_real))

x0+= h

i=1
while(i != n):
    
    dzdx_approx_vector3.append(rungeKutta(x0, dzdx_approx_vector3[i], h, dzdx_func))
    
    z_approx_vector3.append(rungeKutta(x0, z_approx_vector3[i], h, z_func_real))
    
    x0+= h # Update next value of x

    i+=1

print("\n")
print("Vetor X = " + str(x_vector))
print("Vetor Y_Solucao = "  + str(z_approx_vector3))

y_real_sol_vec = []
for j in x_vector:
    y_real_sol_vec.append( m.exp(2)*(m.exp(4) - 1)**(-1) * (m.exp(2*j) - m.exp(-2*j)) + j)


error_vector = [0]
for k in range(n-1):
    error_vector.append( (1 - z_approx_vector3[k+1]/y_real_sol_vec[k+1]) * 100 )

print("\n")
print("Vetor Erro = "  + str(error_vector))
print("\n")
print("Erro Maximo = " + str(max(error_vector)) + "%")

erro_medio = 0
for l in range(n):
    erro_medio += error_vector[k]
erro_medio /= 10

print("Erro Medio = " + str(erro_medio) + "%")


figure, (axis1,axis2,axis3,axis4) = plt.subplots(1,4)

axis1.plot(x_vector,z_approx_vector)
axis1.set_title("Primeiro Chute")

axis2.plot(x_vector,z_approx_vector2)
axis2.set_title("Segundo Chute")

axis3.plot(x_vector,z_approx_vector3)
axis3.set_title("Resultado Aproximado")

axis4.plot(x_vector,y_real_sol_vec)
axis4.set_title("Resultado Real")

plt.show()
