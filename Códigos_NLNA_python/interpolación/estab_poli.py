"""
Estabilidad interpolacion polinomial
"""
from __future__ import division
import numpy as np
from sympy import symbols, prod, lambdify
import matplotlib.pyplot as plt

pi=np.pi

def divdif(x,y):
    n = len(y) 
    d = np.zeros((n , n ))
    d[:,0] = y  # Asigna el vector y a la primera columna de d
    
    for j in range(1, n):
        for i in range(j,n):
            d[i,j] = (d[i-1, j-1] - d[i,j-1]) / (x[i-j] - x[i])
    
    a = np.diag(d) # En la diagonal están los coef f[x0,x1,...,xn]

    z = symbols('z')
    a1 = [z - xi for xi in x]
    a2 = [1] + [prod(a1[:j]) for j in range(1, len(x))]

    # Calcula el polinomio interpolador pn
    pn = sum([ai * ai2 for ai, ai2 in zip(a, a2)])

    # Convierte el polinomio interpolador a una funcion numerica
    pn_f = lambdify(z, pn, modules='numpy')
    return pn_f

def f_stab(x):
    return np.sin(2*pi*x)

# Generación del vector x
x = np.arange(-1, 1.1, 0.1)
x = x[:,np.newaxis]
fs = np.zeros(len(x))

# Llamada a la función 'funcion'
for k in range (len(x)):
    fs[k]=f_stab(x[k])

fs = fs.flatten()[:,np.newaxis]

# Generación de z y w
z = np.random.rand(21, 1)
w = np.dot(np.diag(np.sign(z-0.5).flatten()),9.5e-4*z)

fst = fs + w

# Cálculo de coeficientes y polinomios
pfs = divdif(x.flatten(),fs.flatten())
pfst = divdif(x.flatten(), fst.flatten())

# Evaluación en el rango [-1, 1]
u = np.linspace(-1,1,100)
pfsg = pfs(u)
pfstg = pfst(u)
fr = f_stab(u)

# Gráficos

plt.figure(1)
plt.plot(u, pfsg, '-g', u, fr, ':m')
plt.title('Comparación de pfsg y fr')

plt.figure(2)
plt.plot(u, pfsg, '+b', u, pfstg, '.-r')
plt.title('Comparación de pfsg y pfstg')

plt.show()