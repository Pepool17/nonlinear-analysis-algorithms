"""
Metodo RK 4
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

pi=np.pi
# Righthand side of differential equation
def f(t,y):
    return -y

def yex(t):
    return np.exp(-t)

#Step h
T=10
N=100
h=T/N

# Define discretized time; assumes dt divides nicely into T
t = np.linspace(0.,T,N+1)
#print(t)

# An array to store the solution
eta = np.zeros(len(t))
y=np.zeros(len(t))

# Integrate the differential equation using Euler’s method
eta[0] = 1

for k in range(len(t)-1):
    k1 = f(t[k],eta[k])
    k2 = f(t[k]+0.5*h,eta[k]+0.5*h*k1)
    k3 = f(t[k]+0.5*h,eta[k]+0.5*h*k2)
    k4 = f(t[k+1]+h,eta[k]+h*k3)
    
    eta[k+1] = eta[k]+(1/6)*h*(k1+2*k2+2*k3+k4)

for i in range(0,len(t)):
    y[i] = yex(t[i])
    
error=np.linalg.norm(eta-y, np.inf)
print(error)
# Save the solution
np.savetxt('t.txt', t)
np.savetxt('eta.txt', eta)

plt.plot(t,eta,'+-r',label='Euler explícito')
plt.plot(t,y,'b',label='Exacta')
plt.legend()
plt.show()
# Save figure
#plt.savefig('plot_Euler.pdf')
print(error)

