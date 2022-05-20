# -*- coding: utf-8 -*-
"""
Created on Thu May 19 19:41:55 2022

@author: hp
"""

import matplotlib.pyplot as plt
import numpy as np

def Explicit_pde(boundary, initial, x, y, h, k):
    n = len(x)
    m = len(y)
   # alpha = k/h**2
    U = np.zeros((n,m))
    U[0, :]=boundary[0]
    U[-1, :]=boundary[1]
    U[:, 0]= initial
    print(U.round(3))
    alpha = k/h**2
    for j in range(1,m):
        for i in range(1,n-1):
            U[i,j]=alpha*U[i-1,j-1]+(1-2*alpha)*U[i,j-1]+alpha*U[i+1,j-1]
        
    return U.round(3)
    
h=0.1
k=0.0008
boundary=[0,0] 
x=np.arange(0,2+h, h)
y=np.arange(0,4+k, k)  
initial = 20*abs(np.sin(np.pi*x))
#print(initial)

Z = Explicit_pde(boundary, initial, x, y, h, k)
print("temperature at different grid point",Z)

y1=[0, 0.008, 0.016, 0.04, 0.08, 0.16, 0.4]
for j in range(len(y)): 
        plt.plot(x, Z[:,0])
        plt.plot(x, Z[:,10])
        plt.plot(x, Z[:,20])
        plt.plot(x, Z[:,50])
        plt.plot(x, Z[:,100])
        plt.plot(x, Z[:,200]) 
        plt.plot(x, Z[:,500])
plt.xlabel('distance [m]')
plt.ylabel('Temperature[$\degree$ C]')
plt.legend()
#plt.title('for time step= 10')
plt.legend([f't={value} s' for value in y1])
plt.show()

'''
Initially the temperature was completing one full cycle of |sin(x)| but as
we increase no. of steps we observe at time=0.4 unit it becom half cycle of
that function. Means initially at two point there was source (nearly at x=0.4,1.6 unit) 
but as time passes the temperature at middle start to increase because of
diffusion. Then at t=0.4 unit at the middle near x=1 unit we get maximum temperature.
So as time passes the temperature start to diffuse and from those initial high 
value points and at the middle there is highest tempearture.(little flat curve)
Though the highest value of the temperature decreases as time step increases.
'''
