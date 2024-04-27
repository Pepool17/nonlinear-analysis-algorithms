"""
Código para integración numérica con trapecio: simple y compuesta
"""
from __future__ import division
import numpy as np

pi=np.pi

def f(x):
    return np.exp(x)*np.cos(x)


def trapecio(a,b,f):
    It=0.5*(b-a)*(f(a)+f(b))
    return It

def trapecio_comp(a,b,f,N):
    h=(b-a)/N
    T = np.zeros(N)
    x=np.arange(a,b+h,h)
    for j in range (0,N):
        T[j]=trapecio(x[j],x[j+1],f)
        
    Tc=np.sum(T)
    return Tc


It=trapecio(0,pi,f)
Tc=trapecio_comp(0,pi,f,10000)

error=np.abs(Tc-(-0.5*(np.exp(pi)+1)))
print(It)
print(Tc)
print(error)