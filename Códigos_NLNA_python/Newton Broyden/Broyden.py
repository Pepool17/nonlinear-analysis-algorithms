"""
Código método de Broyden para el problema F(x)=0, con F:Rn ---> Rn, n>1, 
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

def Fun2(xyz):
    x,y,z = xyz
    return [x + y + z - 3,
            x**2 + y**2 + z**2 - 5,
            np.exp(x) + x*y - x*z - 1]

def rosenbrock_function(xy):
    x,y=xy
    return (1 - x)**2 + 100 * (y - x**2)**2

def rosenbrock_gradient(xy):
    x,y=xy
    df_dx = -2 * (1 - x) - 400 * x * (y - x**2)
    df_dy =  200 * (y - x**2)
    return [df_dx, df_dy]

def branin(xy):
    a=1
    b=5.1 / (4 * np.pi**2)
    c=5/pi
    r=6
    s=10
    t=1 / (8 * pi)
    x,y=xy
    return a * (y - b * x**2 + c * x - r)**2 + s * (1 - t) * np.cos(x) + s

def gradiente_branin(xy):
    a=1
    b=5.1 / (4 * pi**2)
    c=5 / np.pi
    r=6
    s=10
    t=1 / (8 * pi)
    x,y=xy
    df_dx = -2 * a * (y - b * x**2 + c * x - r) * (2 * b * x - c) + s * (1 - t) * (-np.sin(x))
    df_dy = 2 * a * (y - b * x**2 + c * x - r) + s * (1 - t) * 0
    return np.array([df_dx, df_dy])


def Broyden(Fun,GFun,x0,Q0):
    nmax = 1000
    x = x0
    x_history = [x0]
    Q=Q0
    err=1
    error=np.zeros(1)
    error[0]=err
    nit=0
    
    Fx = np.array(GFun(x))
    
    while err > epsi and nit < nmax:
        nit=nit+1
        s=LA.solve(Q,-Fx)
        x=x+s
        x_history.append(x.copy())
        err=LA.norm(s,np.inf)
        Fx1=np.array(GFun(x))
        b=Fx1-Fx
        Q=Q+(1/np.dot(s,s))*np.outer(b-Q@s,s)
        Fx = Fx1
        
        error=np.append(error,err)
    
    #x_redondeado = [round(num, 8) for num in x]
    #print(x_redondeado)
    #print(Fun(x_redondeado))
    Fsol=Fun(x)
        # Graficar la evolución de x en cada iteración
    x_history = np.array(x_history).T
    plt.figure(figsize=(10, 6))  # Ajusta el tamaño de la figura según tus preferencias
    for i in range(len(x0)):
        plt.plot(range(nit + 1), x_history[i], label=f'x[{i}]')
    plt.legend()
    plt.xlabel('# iteraciones')
    plt.ylabel('Valor de x')
    plt.title('Evolución de x en cada iteración')
    plt.show()


    plt.figure()
    plt.semilogy(error)
    plt.xlabel('# iteraciones')
    plt.ylabel('||s_k||')
    plt.show()
    return x,nit,err,Fsol

#[x_sol,nit,err,Fsol] = Broyden(rosenbrock_function,rosenbrock_gradient,[0.5,0.5],np.identity(2))
#[x_sol,nit,err,Fsol] = Broyden(Fun,[0,0.5],np.identity(2))
[x_sol,nit,err,Fsol] = Broyden(branin,gradiente_branin,[0.5,0.5],np.identity(2))
print('Broyden')
print('La solucion x* calculada es:', x_sol)
print('El numero de iteraciones es:', nit)
print('El error final es=', err)
print('F(x*)=', Fsol)

#[x_sol,nit,err,Fsol] = Broyden(Fun2,[1.0,1.0,2.0],np.identity(3))
#print('La solución x* calculada es:', x_sol)
#print('El numero de iteraciones es:', nit)
#print('El error final es=', err)
#print('F(x*)=', Fsol)
