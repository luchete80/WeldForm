# https://docs.sympy.org/latest/tutorial/matrices.html
#Plain strain/stress 1 element example
from numpy import *
import numpy.matlib

import array as arr

#X is combines per row to calculate J
p=1.0/1.732050807568877
gauss=[-p,p]
wi=[0.777,0.777]
#Numerated as in Bathe
#X2=matrix([[1,1],[0,1],[0,0],[1,0]])
#Numerated as in deal.ii
X2=matrix([[0,0],[1,0],[0,1],[1,1]])

K=matrix(numpy.matlib.zeros((8, 8)))
B=matrix(numpy.matlib.zeros((3, 8)))
dHxy=matrix(numpy.matlib.zeros((2, 4)))
print (K)
print (X2[0,0])
print (X2[2,1])

ey=200.e9
nu=0.33

ck = ey*(1. - nu) / ((1. + nu)*(1. - 2. * nu))
c=matrix(numpy.matlib.zeros((3, 3)))
c[0,0]=c[1,1]=ck;					
c[0,1]=c[1,0]=ck*nu / (1. - nu)
c[2,2]=ck*(1 - 2 * nu) / (2.*(1. - nu))
print ("C matrix")
print(c)
	
for i in range(2):
    for j in range(2):
        r=gauss[i]
        s=gauss[j]

        #Numerated as in Bathe
        #dHrs=matrix([[(1+s),-(1+s),-(1-s),(1-s)], [(1+r),(1-r),-(1-r),-(1+r)] ])
        #Numerated as in deal.ii
        dHrs=matrix([[-(1-s),(1-s),-(1+s),(1+s)], [-(1-r),-(1+r),(1-r),(1+r)] ])        
        dHrs/=4
        J=dHrs*X2
        dHxy=linalg.inv(J)*dHrs
        for k in range(4):
            B[0,2*k  ]=dHxy[0,k]
            B[1,2*k+1]=dHxy[1,k]
            B[2,2*k  ]=dHxy[1,k]
            B[2,2*k+1]=dHxy[0,k]
        w=0.25
        print (B)
        K+=(B.transpose()*c*B*w)
        print (K)
print (K)
#Boundary conditions
#Numerated as in Bathe
for i in range(8):
    K[4,i] = K[i,4] = 0.0
    K[5,i] = K[i,5] = 0.0
    K[7,i] = K[i,7] = 0.0
	
K[4,4] = K[5,5] = K[7,7] = 1.;
R=[0,0,1000.0,0,0,0,0,0]

# #Numerated as in deal.ii
for i in range(8):
    K[0,i] = K[i,0] = 0.0
    K[1,i] = K[i,1] = 0.0
    K[3,i] = K[i,3] = 0.0
	
K[0,0] = K[1,1] = K[3,3] = 1.;
R=[0,0,0,0,1000.0,0,0,0]

print (K)

U=linalg.solve(K, R)
print ("Results")
print(U)


