# -*- coding: utf-8 -*-
"""
Created on Sun Apr 10 20:22:28 2022

@author: hp
"""
import math

def linConMethod(Xo, m, a, c,
                             n):
 
    # Initialize the seed state
 
    x=[0 for i in range(n)]
    x[0] = Xo
    # Traverse to generate required
    # numbers of random numbers
    for i in range(1, n):
         
        # Follow the linear congruential method
          x[i] = (((x[i - 1] * a) + c) % m)/m
        
    return x


y = linConMethod(1, 572, 16381, 0, 10000)
#print(y)
#y = linConMethod(1, 65, 1021, 0, 1000)                         

def monte_carlo(f, a, b, n):
    #y = linConMethod(1, 572, 16381, 0,10000)
    x=[0 for i in range(n)]
    
    for i in range(n):
        x[i] = y[i]
        result=0.0
        
    for i in range(n):
        result += f(x[i])
        
    result *= (b-a)/n    
    
    return result 
   



# defining the function
f = lambda x: math.sqrt(1-x**2)
print("value of pi by Monte Carlo method", 4*monte_carlo(f, 0, 1, 10000))
print (" {:<8} {:<20} {:<10} ".format('i','pi','N'))

# results are
for i in range (1,101):
     p = 4*monte_carlo(f, 0, 1, 10*i)
     d = {i:[p, 10*i]}
     for  k, v in d.items():
          pi, N= v
          print ("{:<8} {:<20} {:<10} ".format(k, pi, N))   
          
'''
(i) a = 1021, m=65 
value of pi by Monte Carlo method 3.1602257057197076
(ii) a = 16381 and m = 572 
value of pi by Monte Carlo method 3.139458065117981
'''          
         
'''
For MonteCarlo method both (i) a = 1021, m=65 and (ii) a = 16381 and m = 572 
give nearly same result. But (ii) is better as it gives values nearer to 3.14.
I have attached output .txt file for both (i) and (ii).
'''
