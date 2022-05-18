# -*- coding: utf-8 -*-
"""
Created on Thu May 12 17:13:57 2022

@author: hp
"""

import math
import random
# defining the function
f = lambda x: math.e**(-x**2)
  
def monte_carlo(f, a, b, n):
    x=[0 for i in range(n)]
    
    for i in range(n):
        x[i] = random.uniform(a,b)
        result=0.0
        
    for i in range(n):
        result += f(x[i])
        
    result *= (b-a)/n    
    
    return result    

p = monte_carlo(f, 0, 1, 100000)
print("The result of the integration:",p)


'''
The result of the integration: 0.7466359232030313
'''