import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def y(x):
    return np.sin(np.pi*x)
#np.exp(-x)

def f(t,x):
    return 0 
#t*np.cos(np.pi*x)

T = 1  # tiempo total
a = 0  # inicio del dominio
b = 1  # fin del dominio
Nx = 29  # número de puntos en el espacio
Nt = 60 # número de puntos en el tiempo
#alpha = T / Nt 
alpha = 0.001 
h = (b - a) / Nx

x_values = np.linspace(a,b,Nx)
t_values = np.linspace(0, T,Nt )
u = np.zeros((Nt, Nx))
u[0, :] = y(x_values)
u[0,-1] = 0
u[0,0] = 0

for i in range(0,Nt-1):
    for j in range(1,Nx-1):
        u[i+1,j]=alpha/(h*h)*(u[i,j+1]+u[i,j-1]) + u[i,j]*(1-(2*alpha/h**2)) + alpha*f(t_values[i],x_values[j])

# Crear una gráfica tridimensional
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Crear la malla de tiempo y espacio
T_mesh, X_mesh = np.meshgrid(t_values, x_values)

# Graficar la función u en 3D
ax.plot_surface(T_mesh, X_mesh, u.T, cmap='viridis')

# Añadir etiquetas a los ejes
ax.set_xlabel('Tiempo')
ax.set_ylabel('Espacio')
ax.set_zlabel('u')


# Gráficos
for i in range(0, Nt, Nt // 6):  # Graficar cada 1/6 de tiempo
    plt.plot(x_values, u[i, :], label=f'Tiempo: {t_values[i]:.2f}')

plt.xlabel('Espacio (x)')
plt.ylabel('u')
plt.legend()
plt.title('Evolución de u en el tiempo')
plt.show()