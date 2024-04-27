"""
Método de Euler implícito
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

# Funciones dadas
def f(t,y):
    return t*(y**2)#-y

def dfy(t,y):
    return 2*t*y#-1

def yex(t):
    return 2./(2.-t**2) #np.exp(-t)

# Método de Newton
def met_newton_beul(f,dfy,etak,tk,h,M,tol):
    eta = etak
    
    def F(eta):
        return eta - etak - h * f(tk,eta)
    
    def dF(eta):
        return 1 - h * dfy(tk,eta)

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
def Euler_impl(f,dfy,yex,t0,y0,T,n):
    h = T/n
    t = np.linspace(t0,t0+T,n+1)
    eta = np.zeros_like(t)
    
    eta[0] = y0

    for j in range(n):
        eta[j+1] = met_newton_beul(f,dfy,eta[j],t[j+1],h,100,1e-8)

    ye = np.vectorize(yex)(t)  # Uso de np.vectorize para evaluar yex en un array

    error = np.max(np.abs(eta-ye))
    
    w = np.linspace(t0,t0+T,100)
    we = yex(w)

    plt.plot(t, eta, '+r', label='eta=calculada con Euler implícito')
    plt.plot(w, we, 'b', label='y=Solución exacta')
    plt.legend()
    plt.show()

    return t, eta, error, h

[t, eta, error, h] = Euler_impl(f, dfy, yex, 0, 1, 1.3, 1000)
