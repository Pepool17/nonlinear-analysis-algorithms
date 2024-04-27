"""
Código método de la bisección para el problema f(x)=0, con f no lineal
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

pi=np.pi
eps=np.finfo(float).eps

# Función para f(x)=0
def f(x):
    return x-pow(9,-x)#np.cos(2*x)*np.cos(2*x)-x*x#np.log(x)*np.sin(pi*x) 

# Inicializacion: tolerancia, número máximo de iteraciones, intervalo inicial, 
# vector para guardar la historia de la convergencia
epsi=np.sqrt(eps)
nmax=100

a=0
b=1
xk=0.75


fxk=f(xk)

error=np.zeros(1)
    
    
# Lazo principal
err=1
error[0]=err
nit=0
while err > epsi and nit < nmax:
    nit=nit+1
    x=xk-((b-a)/(f(b)-f(a)))*fxk
    err=np.absolute(x-xk)
    xk=x
    fxk=f(xk)
    error=np.append(error,err)

print('La raíz hallada es=', xk)
print('Se alcanzó convergencia en este número de iteraciones=', nit)
print('El error final es=', err)

#Step h
N=1000
h=1/N

# Define discretized time; assumes dt divides nicely into T
t = np.linspace(0,1.5,N+1)

# Evaluación para gráfico de la función en una vecindad de la solución
fun = np.zeros(len(t))
for i in range(0,len(t)):
    fun[i] = f(t[i])

# Posprocesamiento: gráficos

# Gráfico de la función en una vecindad de la solución
plt.figure()
plt.plot(t, fun, color='blue')
plt.plot(xk,fxk, marker="o", color="red")
plt.xlabel('x')
plt.ylabel('f(x)')

# Gráfico de la historia de convergencia
plt.figure()
plt.semilogy(error)
plt.xlabel('# iteraciones')
plt.ylabel('||xk1-xk||')
# Save figure
#plt.savefig('plot_Emod.pdf')