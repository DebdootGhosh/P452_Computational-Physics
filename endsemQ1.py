# -*- coding: utf-8 -*-
"""
Created on Wed May 18 18:08:11 2022

@author: hp
"""

#import random as rd
import math

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
v = linConMethod(0.06, 572, 16381, 0, 1000)

#V=[v for i in range(500)]
#random walk function in 2d so we have 0<=\theta_i <=pi to cover the whole plane . 
#Since each step is of length 1, for step i, Δxi=cos\theta_i and Δy=sin\theta_i
# So \theta_i=2*\pi*r_i where 0<=r_i<=1      
def randwalk(N):
    x=0
    y=0
    X=[]
    Y=[]
    #print('\n-----------RANDOM WALK-----------')
    #print('------------------------------')    
    #print('x   \ty   ')
    #print('------------------------------')
    for i in range(N):
        theta=2*math.pi*v[i]
        #print('%.5f\t%.5f\t'% (x,y) )
        X.append(x)
        Y.append(y)
        x+=math.cos(theta);          
        y+=math.sin(theta); 
    return(X,Y)


#for different random walks of same No. of steps
def randwalk_repeat(N,M):
    X=[]
    Y=[]
    for j in range(M):
        x,y = randwalk(N)
        X.append(x) 
        Y.append(y)
    return X,Y

# calculating average radial distance from the origin for different No. of step 
# calculate rms distance from origin for different No. of steps
#calculate average of  \Delta x, \Delta y

def Measurement(x,y):
    D=0
    R=0
    X=0
    Y=0
    for i in range(500):
        D+=((x[i][len(x[i])-1])**2 + (y[i][len(y[i])-1])**2)**0.5
        R+=((x[i][len(x[i])-1])**2 + (y[i][len(y[i])-1])**2)
        X=X+(x[i][len(x[i])-1])
        Y=Y+(y[i][len(y[i])-1])
    return D/500,(R/500)**0.5,X/500,Y/500

#calling different functions and using input No of steps we will get the answers
X1,Y1=randwalk_repeat(200,500)


#printing results
a,b,c,d=Measurement(X1,Y1)
print("FOR STEPS NO.=200, OVER 500 WALKS ")
print("average rms distance is",b)
print("average displacement in x direction is",c)
print("average displacement in y direction is",d)


'''
FOR STEPS NO.=200, OVER 500 WALKS 
average rms distance is 14.785352665862977
average displacement in x direction is -11.925325353727496
average displacement in y direction is 8.740324288130472

Here we observe the average R_(rms)(N)=14.78 that is nearly (N)^0.5.
'''
