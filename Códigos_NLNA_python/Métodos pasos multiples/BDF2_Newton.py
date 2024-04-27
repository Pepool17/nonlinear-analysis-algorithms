"""
Método BDF2-Newton
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

# Funciones dadas
def f(t,y):
    return -y

def dfy(t,y):
    return -1

def yex(t):
    return np.exp(-t)

# Metodo de Newton
def met_newton_BDF2(f,dfy,etak,tk,etak1,tk1,h,M,tol):
    eta = etak
    
    def F(eta):
        return eta-(4/3)*etak+(1/3)*etak1-(2/3)*h*f(tk,eta)
    
    def dF(eta):
        return 1 - (2/3)*h*dfy(tk,eta)

    err = tol + 1
    nit = 0

    while err > tol and nit < M:
        nit += 1
        diff = F(eta)/dF(eta)
        eta -= diff
        err = abs(diff)

    if nit >= M:
        print('El metodo no converge a la tolerancia deseada y se alcanza el numero maximo de iteraciones.')
    
    return eta

# Metodo de Euler implicito
def BDF2_Newton(f,dfy,yex,t0,y0,T,n):
    h = T/n
    t = np.linspace(t0,t0+T,n+1)
    eta = np.zeros_like(t)
    
    eta[0] = y0
    
    k1 = eta[0]+0.5*h*f(t[0],eta[0])
    eta[1] = eta[0]+h*f(t[0]+0.5*h,k1)

    for j in range(1,len(t)-1):
        eta[j+1] = met_newton_BDF2(f,dfy,eta[j],t[j+1],eta[j-1],t[j-1],h,100,1e-8)

    ye = np.vectorize(yex)(t)  # Uso de np.vectorize para evaluar yex en un array

    error = np.max(np.abs(eta-ye))
    
    w = np.linspace(t0,t0+T,100)
    we = yex(w)

    plt.plot(t, eta, '+r', label='eta=calculada con BDF2')
    plt.plot(w, we, 'b', label='y=Solución exacta')
    plt.legend()
    plt.show()

    return t, eta, error, h

[t,eta,error,h] = BDF2_Newton(f,dfy,yex,0,1,10,100)
print(error)
