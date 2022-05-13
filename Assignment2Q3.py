# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 13:18:30 2022

@author: hp
"""

#import math

def linConMethod(Xo, m, a, c,
                             n):
 
    # Initialize the seed state
 
    x=[0 for i in range(n)]
    x[0] = Xo
    # Traverse to generate required
    # numbers of random numbers
    for i in range(1, n):
         
        # Follow the linear congruential method
          x[i] = -1+2*(((x[i - 1] * a) + c) % m)/m
        
    return x


y = linConMethod(1, 572, 16381, 0, 1000)
#print(y)
#y = linConMethod(1, 65, 1021, 0, 1000)                         

def monte_carlo(f, a, b, n):
    #y = linConMethod(1, 572, 16381, 0,1000)
    x=[0 for i in range(n)]
    
    for i in range(n):
        x[i] = y[i]
        result=0.0
        
    for i in range(n):
        result += f(x[i])
        
    result *= (b-a)/n    
    
    return result 
   



# defining the function for Steinmetz solid
f = lambda x: 4*(1-x**2)
p = monte_carlo(f, -1, 1, 10000)
print("Area of the solid by Monte Carlo", p)

print (" {:<8} {:<20} {:<10} ".format('i','volume','N'))

# results are
for i in range (1,101):
     p = monte_carlo(f, -1, 1, 10*i)
     d = {i:[p, 10*i]}
     for  k, v in d.items():
          pi, N= v
          print ("{:<8} {:<20} {:<10} ".format(k, pi, N))   

'''
(i) a=1021, m=65
Area of the solid by Monte Carlo 5.35123628933519
(ii) a= 16381. m= 572
Area of the solid by Monte Carlo 5.346757586335043
'''
'''
For this question a= 16381. m= 572 gives more accurate result(nearly 5.33)
so I have taken those values. I have attached both output .txt file for 
(i) a=1021, m=65 and (ii) a= 16381. m= 572.
'''
