# https://docs.sympy.org/latest/tutorial/matrices.html
from numpy import *

import array as arr

init_printing(use_unicode=True)
X2=Matrix(zeros(4,2)) 
dHdx=Matrix(zeros(2,4))
J=Matrix(zeros(2,2))

#X is combines per row to calculate J
gauss=[-0.557,0.557]
wi=[0.777,0.777]
X2=[[1,1],[0,1],[0,0],[1,0]]

for i in range(2)
  for j in range(2)
     r=gauss[i]
     s=gauss[j]
#r, s = symbols('r s')

#h=[(1+r)*(1+s)/4,(1-r)*(1+s)/4,(1-r)*(1-s)/4,(1+r)*(1-s)/4]

#Plane Strain 
H=Matrix(zeros(2, 8))
H=([  [h[0],0,h[1],0,h[2],0,h[3],0],[0,h[0],0,h[1],0,h[2],0,h[3]] ])
dHrs=Matrix([ [(1+s)/4,-(1+s)/4,-(1-s)/4,(1-s)/4], [(1+r)/4,(1-r)/4,-(1-r)/4,(1+r)/4] ])
J=dHrs*X2
#B=Matrix(J**(-1)*)
dHxy=linalg.inv(J)*dHrs
for k in range(4)
   B(2*k  ,0)=dHxy(k,0)
   B(2*k+1,1)=dHxy(k,1)
   B(2*k  ,2)=dHxy(k,1)
   B(2*k+1,2)=dHxy(k,0)

K+=(B.transpose()*c*B*w)

#dHdrs is [ [h1r,h2r,h3r,h4r],[h1s,h2s,h3s,h4s])
#dHxy=J-1 dHrs
#B(r,s)=[]
#J  =[[dxr dyr],[dxs dys]]
#J-1=[[drx dsx],[dry dsy]]

#C Matrix 
print(H)
print(dHrs)
print(X2)
print(J)

