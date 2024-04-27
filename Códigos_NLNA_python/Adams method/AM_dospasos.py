"""
Método de Adams-Moulton 2 pasos
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

# Método de Newton
def met_newton_AM_2(f,dfy,etak,etakm1,tk,tk1,tkm1,h,M,tol):
    eta = etak
    
    def F(eta):
        return eta - etak - (h/12)*(5*f(tk1,eta)+8*f(tk,etak)-f(tkm1,etakm1))
    
    def dF(eta):
        return 1 - (h/12)*(5*dfy(tk1,eta))

    err = tol + 1
    nit = 0

    while err > tol and nit < M:
        nit += 1
        diff = F(eta)/dF(eta)
        eta -= diff
        err = abs(diff)

    if nit >= M:
        print('El método no converge a la tolerancia deseada y se alcanzó el número máximo de iteraciones.')
    
    return eta

# Método de Euler implícito
def AM_2(f,dfy,yex,t0,y0,T,n):
    h = T/n
    t = np.linspace(t0,t0+T,n+1)
    eta = np.zeros_like(t)
    
    eta[0] = y0
    # Euler exp. para eta1
    eta[1] = eta[0] + h*f(t[0],eta[0])

    for j in range(1,len(t)-1):
        eta[j+1] = met_newton_AM_2(f,dfy,eta[j],eta[j-1],t[j],t[j+1],t[j-1],h,100,1e-8)

    ye = np.vectorize(yex)(t)  # Uso de np.vectorize para evaluar yex en un array

    error = np.max(np.abs(eta-ye))
    
    w = np.linspace(t0,t0+T,100)
    we = yex(w)

    plt.plot(t, eta, '+r', label='eta=calculada con AM2')
    plt.plot(w, we, 'b', label='y=Solución exacta')
    plt.legend()
    plt.show()

    return t, eta, error, h

[t, eta, error, h] = AM_2(f, dfy, yex, 0, 1, 1.3, 10)
print(error)
