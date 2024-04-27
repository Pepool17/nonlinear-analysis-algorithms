"""
Código método de Newton para el problema f(x)=0, con f no lineal
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

pi=np.pi
eps=np.finfo(float).eps

# Función para f(x)=0
def f(x):
    return np.power(x*x-1,4)*np.log(x)

def df(x):
    return 8*x*np.power(x*x-1,3)*np.log(x)+np.power(x*x-1,4)/x
   

# Inicializacion: tolerancia, número máximo de iteraciones, intervalo inicial, 
# vector para guardar la historia de la convergencia
epsi=np.sqrt(eps)
nmax=100

x=0.8
fx=f(x)
dfx=df(x)

error=np.zeros(1)
r=np.zeros(1)
    
# Lazo principal
err=1
error[0]=err
r[0]=1
m=1

nit=0
while err > epsi and nit < nmax:
    nit=nit+1
    if dfx == 0:
        print('Stop due to vanishing derivative')
    diff=fx/dfx;
    x=x-m*diff
    fx=f(x);
    dfx=df(x);
    rd=err 
    err=np.abs(diff);
    error=np.append(error,err)
    ra=err/rd 
    r=np.append(r,ra) 
    diff=np.absolute(r[nit]-r[nit-1])
    if diff < 1.e-3 and ra > 1.e-2:
        arr=[m,1/(np.absolute(1-ra))]
        m=np.max(arr) 
    else:
        m=m

print('La raíz hallada es=', x)
print('Se alcanzó convergencia en este número de iteraciones=', nit)
print('El error final es=', err)
print('La multiplicidad calculada es=', m)

#Step h
N=1000
h=1/N

# Define discretized time; assumes dt divides nicely into T
t = np.linspace(0.5,1.5,N+1)

# Evaluación para gráfico de la función en una vecindad de la solución
fun = np.zeros(len(t))
for i in range(0,len(t)):
    fun[i] = f(t[i])

# Posprocesamiento: gráficos

# Gráfico de la función en una vecindad de la solución
plt.figure()
plt.plot(t, fun, color='blue')
plt.plot(x,fx, marker="o", color="red")
plt.xlabel('x')
plt.ylabel('f(x)')

# Gráfico de la historia de convergencia
plt.figure()
plt.semilogy(error)
plt.xlabel('# iteraciones')
plt.ylabel('||xk1-xk||')
# Save figure
#plt.savefig('plot_Emod.pdf')