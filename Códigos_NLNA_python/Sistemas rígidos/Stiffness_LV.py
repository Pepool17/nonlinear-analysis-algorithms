"""
Ejemplo rigidez (stiffnes) usando el modelo de Lotka-Volterra

"""
import numpy as np
import matplotlib.pyplot as plt

# Definir las ecuaciones del sistema de Lotka-Volterra
def lotka_volterra(t,y1,y2):
    alpha = 4
    beta = 1
    gamma = 1
    delta = 2
    dy1dt = alpha*y1-beta*y1*y2
    dy2dt = delta*y1*y2-gamma*y2
    return np.array([dy1dt,dy2dt])

# Metodo de Euler
def euler_method(f,t0,y0,T,n):
    h = T/n
    t = np.linspace(t0,t0+T,n+1)
    y = np.zeros((n+1,len(y0)))
    y[0] = y0
    for i in range(n):
        y[i+1] = y[i]+h*f(t[i],*y[i])
    return t, y

# Metodo de Heun
def heun_method(f,t0,y0,T,n):
    h = T/n
    t = np.linspace(t0, t0 + T, n + 1)
    y = np.zeros((n + 1, len(y0)))
    y[0] = y0
    for i in range(n):
        k1 = f(t[i],*y[i])
        k2 = f(t[i]+h,*(y[i]+h*k1))
        y[i+1] = y[i]+0.5*h*(k1+k2)
    return t,y

# Condiciones iniciales y parametros
t0 = 0
T = 15
n = 1000
y0 = np.array([5,5])  # Poblacion inicial de presas y depredadores

# Resolver usando metodo de Euler
[t_euler,y_euler] = euler_method(lotka_volterra,t0,y0,T,n)

# Resolver usando metodo de Heun
[t_heun,y_heun] = heun_method(lotka_volterra,t0,y0,T,n)

# Graficar los resultados
plt.plot(t_euler, y_euler[:,0], label='Presas (Euler)')
plt.plot(t_euler, y_euler[:,1], label='Depredadores (Euler)')
plt.plot(t_heun, y_heun[:,0], label='Presas (Heun)')
plt.plot(t_heun, y_heun[:,1], label='Depredadores (Heun)')
plt.xlabel('Tiempo')
plt.ylabel('Población')
plt.title('Dinamica de presas y depredadores con Euler y Heun')
plt.legend()
plt.grid(True)
plt.show()
