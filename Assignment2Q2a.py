# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 20:38:01 2022

@author: hp
"""

import matplotlib.pyplot as plt
#import math

def linConMethod(Xo, m, a, c, n):
 
    # Initialize the seed state
 
    x=[0 for i in range(n)]
    x[0] = Xo
    # Traverse to generate required
    # numbers of random numbers
    for i in range(1, n):
         
        # Follow the linear congruential method
          x[i] =  (((x[i - 1] * a) + c) % m)/m
        
    return x


'''
v = linConMethod(0.25, 65, 1021, 0, 10000)
k = linConMethod(0.5, 65, 1021, 0, 10000)
'''
#print(v)
#print(y)
#y = linConMethod(1, 65, 1021, 0, 1000) 

v = linConMethod(0.25, 572, 16381 , 0, 10000)
k = linConMethod(0.5, 572, 16381, 0, 10000)

def Throwing(f,N):
    m=0
    for i in range(N):
        x=1*v[i]
        y=1*k[i]
        if f(x,y)<=1:
            m=m+1
    return m/N

# given ellipsoid equation
def circle(x,y):
     return (x)**2+(y)**2
 
#plot the No. of steps Vs area of the circle and compare with analytical value 
def plot_Throwing(f,N):
    x=[]
    y=[]
    
    i=1
    while i<=N:
        x.append(i)
        y.append(4*Throwing(circle, i))
        i=i+1
   
    plt.title("No of steps Vs area")  
    plt.xlabel("No of steps")  
    plt.ylabel("area") 
    plt.plot(x,y,'o-', markersize='2',linewidth='1',label='observed data point')
    plt.legend()
    plt.show()

 
# calling the functions and using the input we will get the answers
print('No. of Steps\tarea ')
print('------------------------------')

for i in range(1,101):
     y=4*Throwing(circle , 100*i)
     print('%.5f\t%.5f\t'% (100*i, y) )

   
plot_Throwing(circle, 10000)  

'''
For throwing method (ii) a = 16381 and m = 572 gives better output(nearly 3.14)
than (i) a = 1021, m=65. So I have taken (ii) for Question no. 3. I have attached 
the output .txt file for both (i) and (ii).
'''