# https://docs.sympy.org/latest/tutorial/matrices.html
from sympy import *

import array as arr

init_printing(use_unicode=True)
X2=Matrix(zeros(4,2)) 
#X is combines per row to calculate J

X2=[[1,1],[0,1],[0,0],[1,0]]

r, s = symbols('r s')
h=[(1+r)*(1+s)/4,(1-r)*(1+s)/4,(1-r)*(1-s)/4,(1+r)*(1-s)/4]

J=Matrix(zeros(2,2))

#Plane Strain 
H=Matrix(zeros(2, 8))
H=([  [h[0],0,h[1],0,h[2],0,h[3],0],[0,h[0],0,h[1],0,h[2],0,h[3]] ])
dHrs=Matrix([ [(1+s)/4,-(1+s)/4,-(1-s)/4,(1-s)/4], [(1+r)/4,(1-r)/4,-(1-r)/4,(1+r)/4] ])
J=dHrs*X2
B=Matrix(J**(-1)*)

#dHdrs is [ [h1x,h2x,h3x,h4x],[h1y,h2y,h3y,h4y])
#dHxy=J-1 dHrs
#B(r,s)=[]
#J  =[[dxr dyr],[dxs dys]]
#J-1=[[drx dsx],[dry dsy]]

#C Matrix (Bathe table 4.3)

expr = sin(r)/r
expr.evalf(subs={r: 3.14})
print (expr.evalf(subs={r: 3.14}))
print (h[1])
print (h[1].evalf(subs={r: 1.0}))
print(H)
print(dHrs)
print(X2)
print(J)

