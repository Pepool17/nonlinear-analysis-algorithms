"""
Código método de Newton para el problema F(x)=0, con F:Rn ---> Rn, n>1, 
no lineal
"""
from __future__ import division
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt

pi=np.pi
eps=np.finfo(float).eps
epsi=np.sqrt(eps)

def Fun(xy):
    x,y = xy
    return [x**2+y**2-1,np.sin(0.5*pi*x)+y**3]

def JFun(xy):
    x,y = xy
    return [[2.*x, 2.*y],
            [0.5*pi*np.cos(0.5*pi*x), 3*y**2]]

# From the exercise:
def Fun2(xyz):
    x,y,z = xyz
    return [x + y + z - 3,
            x**2 + y**2 + z**2 - 5,
            np.exp(x) + x*y - x*z - 1]

def JFun2(xyz):
    x,y,z = xyz
    return [[1, 1, 1],
            [2*x, 2*y, 2*z],
            [np.exp(x) + y - z, x, -x]]


def newton_sist(Fun,x0,JFun):
    nmax = 50

    x = x0
    
    err=1
    error=np.zeros(1)
    error[0]=err
    nit=0

    while err > epsi and nit < nmax:
        nit=nit+1
        Jx = np.array(JFun(x))
        Fx = np.array(Fun(x))
        # 
        # Resolvemos el sistema J(xn)*d_k=-F(xn):
        #
        diff = LA.solve(Jx,-Fx)
        x = x + diff
        
        err=LA.norm(diff)
        error=np.append(error,err)
  
    Fsol=Fun(x)
    plt.figure()
    plt.semilogy(error)
    plt.xlabel('# iteraciones')
    plt.ylabel('||d_k||')

    return x,nit,err,Fsol

#Función 1
[x_sol,nit,err,Fsol] = newton_sist(Fun,[0,0.5],JFun)
print('La solución x* calculada es:', x_sol)
print('El numero de iteraciones es:', nit)
print('El error final es=', err)
print('F(x*)=', Fsol)

# Función 2
[x_sol,nit,err,Fsol] = newton_sist(Fun2,[2.0,1.0,2.0],JFun2)
print('La solución x* calculada es:', x_sol)
print('El numero de iteraciones es:', nit)
print('El error final es=', err)
print('F(x*)=', Fsol)