import numpy as np
import matplotlib.pyplot as plt

# Definir parámetros y funciones
a, b = 0, 1  # Intervalo en x
T = 100  # Tiempo total
N_x, N_t = 100, 9000  # Número de puntos en x y t
h_x = (b - a) / N_x  # Tamaño del paso en x
h_t = T / N_t  # Tamaño del paso en t

def y0(x):
    return  np.exp(-x) 
    #return np.sin(np.pi*x)

def f(t,x):
    return t*np.cos(np.pi*x)


# Inicializar la matriz y con ceros
y = np.zeros((N_t + 1, N_x + 1))

# Establecer la condición inicial
y[0, :] = y0(np.linspace(a, b, N_x + 1))
y[0,0]=0
y[0,-1]=0

# Implementar el método implícito de Euler
for n in range(N_t):
    # Construir la matriz tridiagonal A para el sistema lineal
    A = np.diag(1 + 2 * h_t / h_x**2 * np.ones(N_x - 1)) + np.diag(-h_t / h_x**2 * np.ones(N_x - 2), k=1) + np.diag(-h_t / h_x**2 * np.ones(N_x - 2), k=-1)

    # Definir el lado derecho del sistema lineal
    b_vector = y[n, 1:-1] + h_t * f((n + 1) * h_t, np.linspace(a, b, N_x + 1)[1:-1])

    # Resolver el sistema lineal
    y_n1_interior = np.linalg.solve(A, b_vector)

    # Actualizar la solución para el siguiente paso de tiempo
    y[n + 1, 1:-1] = y_n1_interior

#Crear una malla para plotear
x_mesh, t_mesh = np.linspace(a, b, N_x + 1), np.linspace(0, T, N_t + 1)
x_mesh, t_mesh = np.meshgrid(x_mesh, t_mesh)

# Plotear la solución
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(x_mesh, t_mesh, y, cmap='viridis')
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('y(t, x)')
ax.set_title('Solución de la ecuación diferencial')
plt.show() 
   