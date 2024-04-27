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
    return x-pow(9,-x)
#np.cos(2*x)*np.cos(2*x)-x*x
#np.log(x)*np.sin(pi*x) 

# Inicializacion: tolerancia, número máximo de iteraciones, intervalo inicial, 
# vector para guardar la historia de la convergencia
epsi= np.exp(-10)#np.sqrt(eps)
print(epsi)
nmax=100

a=0#2.5
b=1#3.1

fa=f(a)
fb=f(b)

c=(a+b)/2
fc=f(c)

error=np.zeros(1)

# Control: evitamos iterar si no hay una raíz en [a,b], si la raíz coincide 
# con las fronteras del intervalo o con el punto medio inicial x0=(a+b)/2

if np.sign(fa)==np.sign(fb):
    print('no hay cero de f en [a,b]')
elif np.absolute(fc)<=epsi:
    print('(a+b)/2 es una raiz o esta cerca de una raiz')
    print('raiz=',c)
elif np.absolute(fa)<=epsi: 
    print('a es una raiz o esta cerca de una raiz')
    print('raiz=',a)
elif np.absolute(fb)<=epsi:
    print('b es una raiz o esta cerca de una raiz')
    print('raiz=',b)    
    
    
# Lazo principal
err=b-a
error[0]=err
nit=0
while err > epsi and np.absolute(fc)>epsi and nit < nmax:
    nit=nit+1
    fa=f(a)
    fb=f(b)
    e=0.5*(b-a)
    c=a+e;
    fc=f(c)
    if np.sign(fa)!=np.sign(fc):
        b=c 
    elif np.sign(fc)!=np.sign(fb): 
        a=c
    err=np.absolute(e)
    error=np.append(error,err)

print('La raíz hallada es=', c)
print('Se alcanzó convergencia en este número de iteraciones=', nit)
print('El error final es=', err)

#Step h
N=1000
h=(b-a)/N

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
plt.plot(c,fc, marker="o", color="red")
plt.xlabel('x')
plt.ylabel('f(x)')

# Gráfico de la historia de convergencia
plt.figure()
plt.semilogy(error)
plt.xlabel('# iteraciones')
plt.ylabel('|I_k|')
# Save figure
#plt.savefig('plot_Emod.pdf')
