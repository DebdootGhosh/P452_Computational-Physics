# -*- coding: utf-8 -*-
"""
Created on Thu May 12 21:54:18 2022

@author: hp
"""
import math 
import matplotlib.pyplot as plt
import library as lib

def norm(y,h):
    norm=[0 for i in range(len(y))]
    summ = (lib.dotproduct(y,y)*h)**0.5
    for i in range(len(y)):
          norm[i]=y[i]/summ
    return norm     
# RK-4 method
def rk4(f,x0,y0,v0,xn,h):

    # Calculating no. of steps 
    n = int((xn-x0)/h)
    x=[0 for i in range(n+1)]
    y=[0 for i in range(n+1)]
    # to solve dy/dx = v 
    def v(x,v):
        return v
    '''
    print('\n--------SOLUTION---------------')
    print('--------------------------------')    
    print('x0  \ty0  \tyn ')
    print('--------------------------------')
    '''
    # solving dv/dx = f(x,y,v) and dy/dx=v(x,v) simultaneously
    for i in range(n+1):
        k1 = h * (f(x0, y0,v0))
        k1v = h* (v(x0, v0))
        k2 = h * (f((x0+h/2), (y0+k1/2),v0))
        k2v = h* (v((x0+h/2),v0))
        k3 = h * (f((x0+h/2), (y0+k2/2),v0))
        k3v = h* (v((x0+h/2),v0))
        k4 = h * (f((x0+h), (y0+k3),v0))
        k4v = h* (v((x0+h),v0))
        k= (k1+2*k2+2*k3+k4)/6
        kv = (k1v+2*k2v+2*k3v+k4v)/6
        yn = y0 + kv
        vn = v0 + k
        #print('%.5f\t%.5f\t%.5f'% (x0,y0,yn) )
        #print('-----------------------------')
        x.append(x0)
        y.append(yn)
        y0 = yn
        x0 = x0+h
        v0 = vn
    plt.title("plot of solution for the ground energy state.")  
    plt.ylabel("Wavefunction(Ψ(x))")  
    plt.xlabel("position(x)")
    #plt.plot(x,y,'o', markersize='5',label='given data point')
    plt.plot(x,norm(y,h),'-', markersize='2',linewidth='2',label='solved data')
    plt.legend()
    plt.show()    
    return yn


#shooting method for boundary value problem
def shooting_method(f1,f2,x0,y0,xn,yn,h,g):
     # yn for the first guess
    y = rk4(f2, x0, y0, g, xn, h)
    print(f"Value of y(x=xn) for the above guess {g}=",y)

    # initialising
    lower,upper=0.0,0.0
    # checking if y overshoots or undershoots
    # if y overshoots
    if y > yn and abs(y - yn) > 10**(-4):
        upper = g
        # we get upper bracket of y
        # to find lower bound
        while y > yn:
            g = float(input(f"Guess a value of y\'({x0}) lower than the previous guess\n"))
            y =rk4(f2, x0, y0, g, xn, h)
            print(f"Value of y(x=xn) for the above guess {g}=", y)
         # if yn for the guess is equal to or very near to actual yn   
        if abs(y - yn) < 10E-4:
            y=rk4(f2, x0, y0, g, xn, h)     
            print(f"Value of y(x=xn) for the above guess {g}=", y)
            print("Value of y(x=xn) found, integration successful")
            return 0
        else:                                # if yn of guess is less than actual yn
            lower = g                        # then we have found the lower bracket
            lagrange_interpolation(upper,lower,f1,f2,x0,y0,xn,yn,h)

        # if y undershoots
    elif y < yn and abs(y - yn) > 10E-4:
        lower = g     # got the lower bracket
        # now to find upper bound
        while y < yn:
            g = float(input(f"Guess a value of y\'({x0}) greater than the previous guess\n"))
            y = rk4(f2, x0, y0, g, xn,h)
            print(f"Value of y(x=xn) for the above guess {g}=", y)

        if abs(y - yn) < 10E-4:
            rk4(f2, x0, y0, g, xn, h)      
            print("Value of y(x=xn) found, integration successful")
        else:
            upper = g
            lagrange_interpolation(upper, lower, f1, f2, x0, y0, xn, yn, h)
        # if yn for the guess is equal to or very near to actual yn
    elif abs(y - yn) < 10E-4:           
        y=rk4( f2, x0, y0, g, xn, h)  
        print(f"Value of y(x=xn) for the above guess {g}=", y)
        print("Value of y(x=xn) found, integration successful")



# lagrange interpolation function
def lagrange_interpolation(upper,lower,f1, f2, x0, y0, xn, yn, h):
    # yn for lower bracket
    yl = rk4(f2, x0, y0, lower,xn, h)    
    # yn for upper bracket
    yh = rk4(f2, x0, y0, upper, xn, h)    
    # for next y'(x0)
    g = lower + ((upper - lower) / (yh - yl)) * (yn - yl)
    # yn for the new y'(x0)
    y = rk4(f2, x0, y0, g, xn, h)
    return y
    print("Value of y(x=xn) found, integration successful")
    
    
 
f1=lambda z,y,x:z
f2=lambda z,y,x:-((math.pi)**2)*y
x0,y0=0,0
xn,yn=1,0

h=float(input("Enter step size h: \n"))
print("For h=",h)
print("-----------")
n = int((xn - x0) / h)
g = float(input(f"Guess a value of y\'({x0}) \n"))
shooting_method(f1,f2,x0,y0,xn,yn,h,g)    

'''
Enter step size h: 
0.000001
For h= 1e-06
-----------

Guess a value of y'(0) 
8
Value of y(x=xn) for the above guess 8.0= 1.1739301230999338e-05
Value of y(x=xn) for the above guess 8.0= 1.1739301230999338e-05
Value of y(x=xn) found, integration successful


Enter step size h: 
0.000001
For h= 1e-06
-----------

Guess a value of y'(0) 
33
Value of y(x=xn) for the above guess 33.0= 4.842461719356289e-05
Value of y(x=xn) for the above guess 33.0= 4.842461719356289e-05
Value of y(x=xn) found, integration successful

runfile('C:/Users/hp/new/Assignment3Q2a.py', wdir='C:/Users/hp/new')
Reloaded modules: library


Enter step size h: 
0.000001
For h= 1e-06
-----------

Guess a value of y'(0) 
78
Value of y(x=xn) for the above guess 78.0= 0.00011445818673949584

Guess a value of y'(0) lower than the previous guess
67
Value of y(x=xn) for the above guess 67.0= 9.831664734790418e-05

Guess a value of y'(0) lower than the previous guess
45
Value of y(x=xn) for the above guess 45.0= 6.603356935098236e-05

Guess a value of y'(0) lower than the previous guess
34
Value of y(x=xn) for the above guess 34.0= 4.9892029714596204e-05

Guess a value of y'(0) lower than the previous guess
18
Value of y(x=xn) for the above guess 18.0= 2.6413427676963718e-05

Guess a value of y'(0) lower than the previous guess
8
Value of y(x=xn) for the above guess 8.0= 1.1739301230999338e-05


nter step size h: 
0.000001
For h= 1e-06
-----------

Guess a value of y'(0) 
-8
Value of y(x=xn) for the above guess -8.0= -1.1739301230999338e-05
Value of y(x=xn) for the above guess -8.0= -1.1739301230999338e-05
Value of y(x=xn) found, integration successful


Enter step size h: 
0.00001
For h= 1e-05
-----------

Guess a value of y'(0) 
-88
Value of y(x=xn) for the above guess -88.0= -0.002171457993981092

Guess a value of y'(0) greater than the previous guess
-67
Value of y(x=xn) for the above guess -67.0= -0.0016532691543010513

Guess a value of y'(0) greater than the previous guess
-45
Value of y(x=xn) for the above guess -45.0= -0.0011104046558416192

Guess a value of y'(0) greater than the previous guess
-24
Value of y(x=xn) for the above guess -24.0= -0.0005922158164053554

Guess a value of y'(0) greater than the previous guess
-11
Value of y(x=xn) for the above guess -11.0= -0.0002714322492476365

Guess a value of y'(0) greater than the previous guess
7
Value of y(x=xn) for the above guess 7.0= 0.00017272961309398806
Value of y(x=xn) found, integration successful

'''