# https://docs.sympy.org/latest/tutorial/matrices.html
#Plain strain/stress 1 element example
from numpy import *
import numpy.matlib

import array as arr

#X is combines per row to calculate J
gauss=[-0.557,0.557]
p=1.0/1.732050807568877
wi=[0.777,0.777]
X2=matrix([[1,1],[0,1],[0,0],[1,0]])

K=matrix(numpy.matlib.zeros((8, 8)))
B=matrix(numpy.matlib.zeros((3, 8)))
dHxy=matrix(numpy.matlib.zeros((2, 4)))
print (K)
print (X2[0,0])
print (X2[2,1])

ey=200.e9
nu=0.33

ck = ey*(1. - nu) / ((1. + nu)*(1. - 2 * nu))
c=matrix(numpy.matlib.zeros((3, 3)))
c[0,0]=c[1,1]=ck;					
c[0,1]=c[1,0]=ck*nu / (1. - nu)
c[2,2]=ck*(1 - 2 * nu) / (2.*(1. - nu))
	
	
for i in range(2):
    for j in range(2):
        r=gauss[i]
        s=gauss[j]
#r, s = symbols('r s')

#h=[(1+r)*(1+s)/4,(1-r)*(1+s)/4,(1-r)*(1-s)/4,(1+r)*(1-s)/4]

		#Plane Strain 
		#H=Matrix(zeros(2, 8))
		#H=([  [h[0],0,h[1],0,h[2],0,h[3],0],[0,h[0],0,h[1],0,h[2],0,h[3]] ])
        dHrs=matrix([[(1+s),-(1+s),-(1-s),(1-s)], [(1+r),(1-r),-(1-r),-(1+r)] ])
        J=dHrs*X2
        dHxy=linalg.inv(J)*dHrs
        for k in range(4):
            B[0,2*k  ]=dHxy[0,k]
            B[1,2*k+1]=dHxy[1,k]
            B[2,2*k  ]=dHxy[1,k]
            B[2,2*k+1]=dHxy[0,k]

        K+=(B.transpose()*c*B*w)

R=[0,0,1000.,0,0,0,0,0]
#U=linalg.inv(K)*R

#print(U)

#dHdrs is [ [h1r,h2r,h3r,h4r],[h1s,h2s,h3s,h4s])
#dHxy=J-1 dHrs
#B(r,s)=[]
#J  =[[dxr dyr],[dxs dys]]
#J-1=[[drx dsx],[dry dsy]]

#C Matrix 
#print(H)
print(dHrs)
print(X2)
print(J)

