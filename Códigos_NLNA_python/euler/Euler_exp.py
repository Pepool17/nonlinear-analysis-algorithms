from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

pi=np.pi
# Righthand side of differential equation
def f(x,y):
    #return np.exp(x)*np.sin(2*pi*x) + 2*pi*np.exp(x)*np.cos(2*pi*x)
    return x*(y**2)

def yex(x):
    return 2./(2.-x**2) 
# Solve the differential equation from time 0 to time T
T = 1.4# pi

#Step h
N=1000
h=T/N

# Define discretized time; assumes dt divides nicely into T
t = np.linspace(0.,T,N+1)
#print(t)

# An array to store the solution
eta = np.zeros(len(t))
y=np.zeros(len(t))

# Integrate the differential equation using Euler’s method
eta[0] = 1.
for i in range(1,len(t)):
    eta[i] = eta[i-1] + h*f(t[i-1],eta[i-1])
    #print(eta[i])

for i in range(1,len(t)):
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
print(eta)
