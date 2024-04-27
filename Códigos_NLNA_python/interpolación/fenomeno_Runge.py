"""
Fenómeno Runge
"""
from __future__ import division
import numpy as np
from sympy import symbols, prod, lambdify
import matplotlib.pyplot as plt


def divdif(x,y):
    n = len(y) 
    d = np.zeros((n,n))
    d[:,0] = y  # Asigna el vector y a la primera columna de d
    
    for j in range(1, n):
        for i in range(j,n):
            d[i,j] = (d[i-1, j-1] - d[i,j-1])/(x[i-j]-x[i])
    
    a = np.diag(d) # En la diagonal están los coef f[x0,x1,...,xn]

    z = symbols('z')
    a1 = [z - xi for xi in x]
    a2 = [1] + [prod(a1[:j]) for j in range(1, len(x))]

    # Calcula el polinomio interpolador pn
    pn = sum([ai * ai2 for ai, ai2 in zip(a, a2)])

    # Convierte el polinomio interpolador a una funcion numerica
    pn_f = lambdify(z, pn, modules='numpy')
    return pn_f

def f_Runge(x):
    return 1/(1+x**2)

x5 = np.linspace(-5,5,5)
x10 = np.linspace(-5,5,10)
x15 = np.linspace(-5,5,15)

f5 = np.zeros(len(x5))
f10 = np.zeros(len(x10))
f15 = np.zeros(len(x15))


for k in range (5):
    f5[k]=f_Runge(x5[k])
    
for k in range (10):
    f10[k]=f_Runge(x10[k])

for k in range (15):
    f15[k]=f_Runge(x15[k])    

p5=divdif(x5,f5)
p10=divdif(x10,f10)
p15=divdif(x15,f15)

xg=np.linspace(-5,5,100)
f_R=np.zeros(len(xg))
pol5=np.zeros(len(xg))
pol10=np.zeros(len(xg))
pol15=np.zeros(len(xg))
    
for k in range (len(xg)):
    f_R[k]=f_Runge(xg[k])
    
for k in range (len(xg)):
    pol5[k]=p5(xg[k])
    
for k in range (len(xg)):
    pol10[k]=p10(xg[k])

for k in range (len(xg)):
    pol15[k]=p15(xg[k])
    
plt.plot(xg, f_R, 'b', xg,pol5,'r', xg,pol10,'g')#,xg,pol15,'k')  