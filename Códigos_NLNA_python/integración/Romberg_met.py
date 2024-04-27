"""
Código para integración numérica con método de Romberg
"""
from __future__ import division
import numpy as np

pi=np.pi

def f(x):
    return np.exp(x)*np.cos(x)

def trapecio_comp(a,b,f,N):
    h=(b-a)/N
    T = np.zeros(N)
    x=np.arange(a,b+h,h)
    for j in range (0,N):
        T[j]=0.5*(x[j+1]-x[j])*(f(x[j])+f(x[j+1]))
        
    Tc=np.sum(T)
    return Tc

def Romberg_met(a,b,f,n):
    A = np.zeros((n+1,n+1))
    
    for i in range(n+1):
        A[i,0] = trapecio_comp(a,b,f,2**i)
    
    for j in range(1,n+1):
        for i in range(j,n+1):
            A[i,j] = (4**j*A[i,j-1]-A[i-1,j-1])/(4**j-1)
    
    return A[n,n]

n=8
error=np.zeros(n)

for l in range (n):
    error[l]=np.abs(Romberg_met(0,pi,f,l)-(-0.5*(np.exp(pi)+1)))
    
result = Romberg_met(0,pi,f,n)
print(f"Aproximación de la integral: {result}")
print(error)