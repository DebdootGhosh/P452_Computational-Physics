# -*- coding: utf-8 -*-
"""
Created on Fri May 13 17:49:49 2022

@author: hp
"""

import matplotlib.pyplot as plt
import numpy as np
import library as lib

def diff(U,V):
    X = lib.matrix_subtraction(U,V)
    norm=0
    for i in range(len(X)):
        for j in range(len(X[0])):
            norm = norm+X[i][j]
    return abs(norm)/(len(X))**2
        
def Laplace(Vtop, Vbottom, Vleft, Vright, xmax, ymax, dimnx, dimny):
    
    X,Y = np.meshgrid(np.arange(0, dimnx),np.arange(0,dimny))
    V = [[0 for i in range(dimnx)]for j in range(dimny)]
    Vnew = [[0 for i in range(dimnx)]for j in range(dimny)]   
    for j in range(dimny):
        V[dimny-1][j]=Vtop
    for j in range(dimny):    
         V[0][j]=Vbottom
    for i in range(dimnx):
        V[i][dimnx-1]=Vright
    for i in range(dimnx):
         V[i][0]=Vleft
    #print(diff(V,Vnew))   
    
    k=0
    while diff(V,Vnew) > 0.000001:
        Vnew = np.copy(V)
        for i in range(1,dimnx-1):
            for j in range(1,dimny-1):
                 V[i][j]=0.25*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1])
        #print(k)
        k=k+1
        
    return V
    
    
Z = Laplace(0, 1, 0, 0, 1, 1, 10,10)
print("Voltage at different grid point",lib.print_matrix(Z)) 
x=np.linspace(0,1,len(Z))
x,y=np.meshgrid(x,x)
ax = plt.axes(projection='3d')
Z=np.array(Z)
ax.plot_surface(x, y, Z, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('V')
ax.set_title('3D surface Plot')
plt.show();

colourMap = plt.cm.jet
colorinterpolation = 50
x=np.linspace(0,1,len(Z))
x,y=np.meshgrid(x,x)
plt.title("Contour 2d plot for Voltage")  
plt.contourf(x,y,Z, colorinterpolation, cmap=colourMap)       
plt.colorbar()
plt.show()    