"""
Codigo en DF para resolver el problema -u''(x)+ q(x)*u(x) = f(x)
con condiciones Dirichlet homogeneas. 
@author: sergiogonzalez
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA

pi=np.pi
#
# Funciones f, q y sol exacta
#
def f(x):
    return (1+pi**2)*np.sin(pi*x)

def q(x):
    return 1

def sol(x):
    return np.sin(pi*x)

#
# rutina para armar la matriz de DF
#
def build_matrix_df(n):
    size = n-1
    matrix = np.zeros((size, size), dtype=float)
    for i in range(0,n-1):
        matrix[i,i] = 2. + h**2*qx[i]
    for i in range(1,n-1):
        matrix[i,i-1] = -1.
        matrix[i-1,i] = -1.
    return matrix

#
# tamaño de malla, paso de discretizacion
#
n=5
L=1
h = L/n;

#
# inicializacion vectores
#
x = np.linspace(0.,L,n+1)
xg = np.linspace(0.,L,100)

fx = np.zeros(len(x)-2)
qx = np.zeros(len(x)-2)
u = np.zeros(len(x))
ug = np.zeros(len(xg))
ue = np.zeros(len(x))
err = np.zeros(len(x))

#
# construccion vectores q y f y matriz Ah
#

for i in range(0,len(x)-2):
    fx[i] = f(x[i+1])

for i in range(0,len(x)-2):
    qx[i] = q(x[i+1])

Ah=(1/h**2)*build_matrix_df(n)
print(Ah)
print(fx)
#
# solucion del sistema Ah y = f
#
y = LA.solve(Ah,fx)

#
# solución con C.F. y solucion exacta para grafico
# 
for i in range(1,len(x)-1):
    u[i]=y[i-1]
    
for i in range(0,len(xg)):
    ug[i]=sol(xg[i])
    
#
# error en norma inf
#
for i in range(0,len(x)):
    ue[i]=sol(x[i])
    
for i in range(0,len(x)):
    err[i]=u[i]-ue[i]

error=LA.norm(err,np.inf)
print(error)    

np.savetxt('x.txt', x)
np.savetxt('u.txt', u)

# Create a new figure
plt.figure()
# Plot solution
plt.plot(x,u, color='blue', marker='.')
plt.plot(xg,ug,color='red')
# Axis labels
plt.xlabel('x')
plt.ylabel('u(x)')
# Save figure
plt.savefig('plot_u_SL.pdf')    
#print(u)