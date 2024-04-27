"""
Codigo para calcular las diferencias divididas de Newton
"""
from __future__ import division
import numpy as np
from sympy import symbols, prod

def divdif(x,y):
    n = len(y) 
    d = np.zeros((n,n))
    d[:,0] = y  # Asigna el vector y a la primera columna de d
    
    for j in range(1, n):
        for i in range(j,n):
            d[i,j] = (d[i-1, j-1] - d[i,j-1]) / (x[i-j] - x[i])
    
    a = np.diag(d) # En la diagonal están los coef f[x0,x1,...,xn]
    return a,d

def poli_sym(a, x):
    z = symbols('z')
    a1 = [z - xi for xi in x]
    a2 = [1] + [prod(a1[:j]) for j in range(1, len(x))]

    # Calcula el polinomio interpolador pn
    pn = sum([ai * ai2 for ai, ai2 in zip(a,a2)])

#    pn = sum([ai * ai2 for ai, ai2 in zip(a, a2)])
    
    return pn

# Ejemplo de uso:
x = np.array([0, 0.5, 1])
y = np.array([-1, 1/6, 8/9])  # Ahora y es un vector
[a,d] = divdif(x,y)

pn=poli_sym(a,x)
print('a:', a)
print('pn(z)=',pn)