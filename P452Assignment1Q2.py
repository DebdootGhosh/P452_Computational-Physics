# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 10:52:15 2022

@author: hp
"""

import library as lib
import matplotlib.pyplot as plt
b =[-5/3, 2/3, 3, -4/3, -1/3, 5/3]
A = open('matrix2.txt','r')
m = lib.readwritematrix(A)
print('The given A matrix:',lib.print_matrix(m))
print('print the b vector:', b)

B = open('new_matrix2.txt','r')
j = lib.readwritematrix(B)
print('The given B matrix(for conjugate gradient):',lib.print_matrix(j))


print('Linear System solution by jacobi iteration method:',lib.jacobi_iteration(m, b))
print('Linear System solution by LU decomposition method:',lib.LUsolution(m, b))
#print(lib.seidel(m, b))
#print(lib.conjugate_grad(m, b))
def inverse_jacobi_iteration(A):
    n = len(A)
    x=[0 for i in range(n)]
    inv = []
    for i in range(n):
        x[i] = 1
        y = lib.jacobi_iteration(A, x)
        inv.append(y)
        x[i]= 0
    return lib.transpose(inv)

def inverse_seidel(A):
    n = len(A)
    x=[0 for i in range(n)]
    inv = []
    for i in range(n):
        x[i] = 1
        y = lib.seidel(A, x)
        inv.append(y)
        x[i]= 0
    return lib.transpose(inv)

def inverse_conjugate_grad(A):
    n = len(A)
    x=[0 for i in range(n)]
    inv = []
    for i in range(n):
        x[i] = 1
        y = lib.conjugate_grad(A, x)
        inv.append(y)
        x[i]= 0
    return lib.transpose(inv)

print('The inverse of the matrix A by Jacobi iteration method is:', lib.print_matrix(inverse_jacobi_iteration(m)))
print('The inverse of the matrix A by Gauss Seidel method is:', lib.print_matrix(inverse_seidel(m)))
print('The inverse of the matrix A by Conjugate Gradient method is:', lib.print_matrix(inverse_conjugate_grad(j)))

h = [ 0.65054339651863, 0.4145556547956528,0.2339389456697146, 0.09606181004682304,0.06007867121112899,0.05822342789540097
,0.031196895867116547,0.024888351438532956,0.02518024490001202,0.012289445531384474
,0.0059074522228674264,0.004360787700463421,0.0003449459394575864,0.0003763864485171549
,0.000804293236896398,0.00011007731629219524,2.6716009611143793e-05
,2.699688779939936e-07,5.295947517383195e-06,4.7616771641742205e-07
,1.133670939856354e-06,2.7356947572437797e-06,9.765406779363904e-07]

g =[1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]

plt.title("residue Vs iteration no. plot for Gauss Seidel method (Q2)")  
plt.xlabel("iteration no.")  
plt.ylabel("residue")
plt.plot(g,h,'o', markersize='5',label='residue data point')
plt.plot(g,h,'-', markersize='2',linewidth='2',label='curve')
plt.legend()
plt.show()
    
m = [0.5964541782123884,0.42221689266904183,0.2989023706250272
,0.13283363848936203,0.07597468726708412,0.05197671879505428,0.045971406728955305,0.03090392209474427
,0.019476940432410897,0.018171942631368593,0.01880657742236788
,0.0026666505728952248,0.0029378787833725228,0.0008879360994424558
,0.004133645220788585,0.0026543302956184987,0.0010854105684635616
,0.0004292769359014614,0.00010111639384513282,2.9407532721928106e-05
,2.297060323299559e-05,8.552124443849345e-06,7.910884290677708e-06
,1.264686737984215e-05,1.0097254724920685e-06,4.106440390427084e-06
,6.419470597250385e-08,3.4781494632743702e-06]

n =[1, 2, 3, 4, 5, 6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28]

plt.title("residue Vs iteration no. plot for Jacobi iterative method (Q2)")  
plt.xlabel("iteration no.")  
plt.ylabel("residue")
plt.plot(n,m,'o', markersize='5',label='residue data point')
plt.plot(n,m,'-', markersize='2',linewidth='2',label='curve')
plt.legend()
plt.show()


v = [8.553537713032693,3.0730214765396155,1.847777777777777
,0.4835303119852634,1.5549958370982542e-30]
k = [1,2,3,4,5]

plt.title("residue Vs iteration no. plot for Conjugate Gradient method (Q2)")  
plt.xlabel("iteration no.")  
plt.ylabel("residue")
plt.plot(k,v,'o', markersize='5',label='residue data point')
plt.plot(k,v,'-', markersize='2',linewidth='2',label='curve')
plt.legend()
plt.show()

'''
 I have taken the given matrix for the solution of the
 linear system as well as for the matrix inversion by Gauss-Seidel 
 and Jacobi method. But for conjugate gradient method I have taken 
 a symmetric matrix by changing the given matrix little bit. 
 For Gauss- sedidel and Jacobi , the inverse of the matrix is 
 almost same and their residue vs iteration no. plots are nearly 
 same.
'''
