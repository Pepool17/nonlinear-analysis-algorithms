import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def y(x):
    # return np.sin(np.pi*x)
    return np.exp(-x)

def f(t, x):
    #return 0
    return t * np.cos(np.pi * x)

T = 100  # tiempo total
a = 0  # inicio del dominio
b = 1  # fin del dominio
Nx = 500  # número de puntos en el espacio
Nt = 100 # número de puntos en el tiempo
alpha = 0.001  # Ajusta este valor según la estabilidad
h = (b - a) / Nx

x_values = np.linspace(a, b, Nx)
t_values = np.linspace(0, T, Nt)
u = np.zeros((Nt, Nx))
u[0, :] = y(x_values)

for i in range(0, Nt-1):
    for j in range(1, Nx-1):
        u[i+1, j] = (alpha * (u[i, j+1] + u[i, j-1]) + (2 - 2*alpha) * u[i, j]
                     + alpha * h**2 * f(t_values[i], x_values[j])) / (1 + alpha * h**2)

# Graficar resultados
X, T = np.meshgrid(x_values, t_values)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, T, u, cmap='viridis')
ax.set_xlabel('Espacio (x)')
ax.set_ylabel('Tiempo (t)')
ax.set_zlabel('u(x, t)')
ax.set_title('Evolución de la ecuación de difusión')
plt.show()