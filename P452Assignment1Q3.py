# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 15:32:07 2022

@author: hp
"""

import library as lib
import math
import matplotlib.pyplot as plt

def Af(x):
    m = len(x)
    n = int(math.sqrt(m))
    m2 =  0.04
    a=[]
    for k in range(m):
         i = int(k/n)
         j = k % n
    # Periodic boundary condition
         if i>=n-1: x1= n*(i+1-n)+j
         else: x1 = n*(i+1)+j
    
         if j>=n-1: y1= n*i+(j+1-n)
         else: y1 = n*i+(j+1)
    
         if i<=0: x2= n*(i-1+n)+j
         else: x2 = n*(i-1)+j
    
         if j<=0: y2 = n*i+(j-1+n)
         else: y2 = n*i+(j-1)
         
        
         A1 = 0.5*(x[x1]+x[x2]-2*x[k]+x[y1]+x[y2]-2*x[k])+m2*x[k]
         
         a.append(A1)
         
    return a    
   

# conjugate gradient method
def conjugate_grad( b, x=None):
    """
    Description
    -----------
    Solve a linear equation Ax = b with conjugate gradient method.
    Parameters
    """
    n = len(b)
    if not x:
        x = [1 for _ in range(len(b))]
    r =lib.vector_subtraction(Af(x),b)
    p = [-r[i] for i in range(n)]
    r_k_norm = (lib.dotproduct(r, r))
    for i in range(2*n):
        Ap = Af(p)
        alpha = r_k_norm / lib.dotproduct(p, Ap)
        x = [x[i]+(alpha * p[i]) for i in range(n)]
        r = [r[i]+(alpha * Ap[i]) for i in range(n)]
        r_kplus1_norm = (lib.dotproduct(r, r))
        beta = r_kplus1_norm / r_k_norm
        residue = r_kplus1_norm
        h = []
        g =[]
        h.append(residue)
        g.append(i)
        #print(h)
        r_k_norm = r_kplus1_norm
        if r_kplus1_norm < 1e-6:
            break
        p = [(beta * p[i]) - r[i] for i in range(n)]   
    print('Iteration No.', i)
    return x    


def inverse_conjugate_grad(n):
    x=[0 for i in range(n)]
    inv = []
    for i in range(n):
        x[i] = 1
        y = conjugate_grad(x)
        inv.append(y)
        x[i]= 0
    return lib.transpose(inv)

print('The inverse of the matrix A by Conjugate Gradient method is:', lib.print_matrix(inverse_conjugate_grad(400)))

h = [77.41567947391056,6.928686949301374,2.294906092093678
,1.0314216248644092,0.5733846852616792,0.3522546202490805
,0.24041823081417787,0.17464575490223197,0.13701944744421637
,0.11770020896849002,0.12010875517534328,0.14088253676308204
,0.13975683380331833,0.09564077746919888,0.07119041828672154
,0.05920469024324161,0.0681839871778778,0.1010789654887527
,0.11010197200794324,0.10172786189225015,0.08766880303139757
,0.07100059069756066,0.05663928508907445,0.03840607465073095
,0.022755558460305803,0.011095991999905418,0.004825180083207987
,0.001235449487927272,0.00021638451754166165,4.285518497874467e-05
,1.546288748680275e-05,5.29333726799645e-06,1.9069236296732834e-06
,5.316252122474021e-07]

g =[1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,
    22,23,24,25,26,27,28,29,30,31,32,33,34]

plt.title("residue Vs iteration no. plot for Conjugate Gradient method (Q3)")  
plt.xlabel("iteration no.")  
plt.ylabel("residue")
plt.plot(g,h,'o', markersize='5',label='residue data point')
plt.plot(g,h,'-', markersize='2',linewidth='2',label='curve')
plt.legend()
plt.show()

