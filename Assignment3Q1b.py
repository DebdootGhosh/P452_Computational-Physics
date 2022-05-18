# -*- coding: utf-8 -*-
"""
Created on Thu May 12 18:14:01 2022

@author: hp
"""

import math
import random
# defining the function
f = lambda x: math.e**(-x**2)

def pdf(x,alpha):
    return alpha*math.e**(-x)

alpha = math.e/(math.e-1)

def monte_carlo(f, a, b, n):
    x=[0 for i in range(n)]
    y=[0 for i in range(n)]
    for i in range(n):
        x[i] = random.uniform(a,b)
        y[i] = -math.log(1-x[i]/alpha)
        result=0.0
        
    for i in range(n):
        result += f(y[i])/pdf(y[i],alpha)
        
    result *= (b-a)/n    
    
    return result    

p = monte_carlo(f, 0, 1, 100000)
print("The result of the integration:",p)

'''
The result of the integration: 0.7468239025674249
'''

'''
So this importance sampling method sampling method provides better result 
than simple Monte Carlo method. (near to calculated value of the integration)
'''