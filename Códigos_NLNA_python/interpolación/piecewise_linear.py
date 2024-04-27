"""
Código para interpolación lineal por partes
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

pi=np.pi
eps=np.finfo(float).eps
epsi=np.sqrt(eps)

def f(x):
    return 1/(1+x**2)

a=-5
b=5
N=10
h=(b-a)/N
x=np.zeros(N+2)
y=np.zeros(N+2)

for k in range (0,N+2):
    x[k]=a+k*h

for i in range (0,N+2):
    y[i] = f(x[i])

l = len(x)
m = (y[1:l]-y[0:(l-1)])/(x[1:l]-x[0:(l-1)])
c = (y[0:(l-1)]*x[1:l]-y[1:l]*x[0:(l-1)])/(x[1:l]-x[0:(l-1)])

z = x[0:l-1]
p = m*z+c

print(z,m,p)

w = np.linspace(a,b,num=100)
u=np.zeros(len(w))
for i in range (0,len(w)):
    u[i] = f(w[i])

plt.plot(z, p, 'o-b', w, u, 'r')

