# -*- coding: utf-8 -*-
"""
Created on Sun Apr  3 10:29:41 2022

@author: hp
"""

import library as lib
b =[19, 2, 13, -7, -9, 2]
A = open('matrix1.txt','r')
m = lib.readwritematrix(A)
print('The given A matrix:',lib.print_matrix(m))
print('print the b vector:', b)


print('Linear System solution by Gauss-Jordon method:',lib.gaussjordan(m, b))
print('Linear System solution by LU decomposition method:',lib.LUsolution(m, b))

'''
The given A matrix:
[1.0, -1.0, 4.0, 0.0, 2.0, 9.0]
[0.0, 5.0, -2.0, 7.0, 8.0, 4.0]
[1.0, 0.0, 5.0, 7.0, 3.0, -2.0]
[6.0, -1.0, 2.0, 3.0, 0.0, 8.0]
[-4.0, 2.0, 0.0, 5.0, -5.0, 3.0]
[0.0, 7.0, -1.0, 5.0, 4.0, -2.0]

print the b vector: [19, 2, 13, -7, -9, 2]

Linear System solution by Gauss-Jordon method: [-1.7618170439978567, 0.8962280338740136, 4.051931404116157, -1.6171308025395428, 2.041913538501914, 0.15183248715593495]

Linear System solution by LU decomposition method: [-1.7618170439978567, 0.8962280338740136, 4.051931404116157, -1.6171308025395428, 2.041913538501914, 0.15183248715593495]
'''