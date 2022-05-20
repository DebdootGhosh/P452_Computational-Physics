# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 17:23:06 2022

@author: hp
"""

import library as lib
import matplotlib.pyplot as plt
import numpy as np

def phi0(x):
    return 1

def phi1(x):
    return x

def phi2(x):
    return (1/2)*((3*(x)**2)-1)

def phi3(x):
    return (1/2)*((5*(x)**3)-3*x)
def phi4(x):
    return (1/8)*((35*x**4)-(30*x**2)+3)
def phi5(x):
    return (1/8)*((63*x**5)-(70*x**3)+15*x)#phi5(X[i])
def phi6(x):
    return (1/16)*((231*x**6)-(315*x**4)+(105*x**2)-5)
def phi(X):
    a=[]
    for i in range(len(X)):
        a1=[phi0(X[i]),phi1(X[i]),phi2(X[i]),phi3(X[i]),phi4(X[i]),phi5(X[i]),phi6(X[i])]
        a.append(a1)
    return lib.transpose(a)   
        
    
def foo(x,c):
    z=[]
    
    for i in range(len(x)):
        f=0
        for j in range(len(c)):
            f = f+c[j][0]*phi(x)[j][i]
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
                a1 = power(power(phi(x)[j], 1, phi(x)[i], 1), 1, w(sig), 1)
                a.append(sumv(a1))
                #print(a)
        b1=power(power(phi(x)[i], 1, y, 1), 1, w(sig), 1)
        b.append(sumv(b1))
        #print(b)
        A.append(a)
    C=lib.inverse(A)
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



A = open('esem4fit.txt','r')
m = lib.readwritematrix(A)
n = lib.transpose(m)
x = n[0]
y = n[1]
sig = [1 for i in range(len(x))]


k,c = polyfit(x, y, sig, 4)
#print(len(x))
#print(len(y))

print('Array of coefficient is:',k)
n = lib.transpose(k)
z = foo(x,n)
#print('A matrix is:',c)
#print('Y data from fiiting:',z)
plt.title("plot of fitting the data set with Legendre basis(Q2)")  
plt.xlabel("Y")  
plt.ylabel("X")
plt.plot(x,y,'o', markersize='5',label='given data point')
plt.plot(x,z,'-', markersize='2',linewidth='2',label='fitted curve for n=4')
plt.legend()
plt.show()
f = np.linalg.cond(c)
print('condition number',f)

'''
for order n=4 of Legendre polynomial 
Array of coefficient is: [0.0696577968718434, 0.0036240203429267976, -0.012082580199522458, 0.011426217647052338, 0.11049235140900118]
condition number 7.48446235910864

for order n=5 of Legendre polynomial 
Array of coefficient is: [0.06965779687198835, 0.004301685837864658, -0.012082580199514597, 0.013083743602879054, 0.11049235140900898, -0.006726972223322121]
condition number 9.389111166327133

for order n=6 of Legendre polynomial
Array of coefficient is: [0.07003196671661188, 0.004301685837864297, -0.010166710608983586, 0.013083743602879169, 0.11411855049268517, -0.006726972223322203, -0.012384559712843649]
condition number 11.501156721269945


So order n=4 of Legendre polynomial does the fitting most efficiently as the 
condition number is least for it among n=4, 5, 6. The plots are 
attached in classroom.
'''