"""
Created on Thu Jan 18 14:26:56 2024
Método del punto medio para y'(t)=f(t,y(t)), y(t0)=y0
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

pi=np.pi
# Righthand side of differential equation
def f(x,y):
    #return np.exp(x)*np.sin(2*pi*x) + 2*pi*np.exp(x)*np.cos(2*pi*x)
    return -y#x*(y**2)

def yex(x):
    return np.exp(-x)#2./(2.-x**2) 
# Solve the differential equation from time 0 to time T
T = 1.4# pi

#Step h
N=100
h=T/N

# Define discretized time; assumes dt divides nicely into T
t = np.linspace(0.,T,N+1)
#print(t)

# An array to store the solution
eta = np.zeros(len(t))
y=np.zeros(len(t))

# Condicion inicial
eta[0] = 1.

# Euler exp. para eta1
eta[1] = eta[0] + h*f(t[0],eta[0])

for k in range(1,len(t)-1):
    eta[k+1]=eta[k-1]+2*h*f(t[k],eta[k])

for i in range(0,len(t)):
    y[i] = yex(t[i])
    
error=np.linalg.norm(eta-y, np.inf)
print(error)
# Save the solution
np.savetxt('t.txt', t)
np.savetxt('eta.txt', eta)

plt.plot(t,eta,'+-r',label='Punto medio')
plt.plot(t,y,'b',label='Exacta')
plt.legend()
plt.show()
# Save figure
#plt.savefig('plot_Euler.pdf')

