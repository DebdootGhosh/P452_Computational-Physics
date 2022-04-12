# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 15:22:11 2022

@author: hp
"""
import library as lib
import matplotlib.pyplot as plt
#import math 
import numpy as np

def foo(x,c):
    z=[]
    
    for i in range(len(x)):
        f=0
        for j in range(len(c)):
            f = f+c[j][0]*pow(x[i],j)
        z.append(f)
    return z

def sumv(x):
    Sum=0
    for i in range(len(x)):
        Sum = Sum+x[i]
    return Sum
    
# FOR x^a y^b type value generation
def power(x,a,y,b):
    z=[]
    for i in range(len(x)):
        z.append(pow(x[i],a)*pow(y[i],b))
    return z

def w(sig):
    w=[]
    for i in range(len(sig)):
        w.append((1/sig[i])**2)
    return w

def polyfit(x,y,sig,n):
    b=[]
    A=[]
    # n = order of fit
    for i in range(n+1):
        a=[]
        for j in range(n+1):
            if i==0 and j==0:
                #sumv(w(sig))
                a.append(sumv(w(sig)))
            else:
                a1 = power(power(x, i+j, y, 0), 1, w(sig), 1)
                a.append(sumv(a1))
        b1=power(power(x, i, y, 1), 1, w(sig), 1)
        b.append(sumv(b1))
        A.append(a)
    C=lib.inverse(A)
    print(C)
    X=lib.multiply(C, b)
    return X,C

def Error(x,y,sig,n):
    b=[]
    k,c = polyfit(x, y, sig, n)
    for i in range(len(c)):
        b.append(c[i][i])
        
    return b

def Chi(x,y,sig,n):
    k,c = polyfit(x, y, sig, n)    
    z = foo(x,k) 
    z1 = lib.vector_subtraction(y, z) 
    ch = power(z1,2,w(sig),1)
    
    return sumv(ch)/(len(x)-(n+1))

x = [0.00,0.05 ,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65 	
,0.70,0.75,0.80,0.85,0.90,0.95,1.00]

y = [0.486,0.866,0.944,1.144,1.103,1.202,1.166,1.191,1.124,1.095,1.122,1.102
 	,1.099,1.017,1.111,1.117,1.152,1.265,1.380,1.575,1.857]
sig = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
k,c = polyfit(x, y, sig, 3)
print('Array of coefficient is:',k)
n = lib.transpose(k)
z = foo(x,n)
print('A matrix is:',c)
print('Y data from fiiting:',z)
plt.title("plot of cubic least square fit(Q1)")  
plt.xlabel("Y")  
plt.ylabel("X")
plt.plot(x,y,'o', markersize='5',label='given data point')
plt.plot(x,z,'-', markersize='2',linewidth='2',label='fitted curve')
plt.legend()
plt.show()
f = np.linalg.cond(c)
print('condition number',f)

'''
The condition number of a matrix can indicate how well-conditioned 
the matrix is - if it is very large, the matrix is ill-conditioned, and 
the inversion is prone to errors. But for this original basis we get condition
no of the matrix is very larger than modified basis so the modified basis 
functions better fit the given data.
'''