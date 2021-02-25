from Domain import *
from Vec3_t import *
from Particle import *

dom= Domain()
x=Vec3_t (0,0,0)
print(x)
print ("Hola")

dom.Dimension	= 2
dom.Nproc	    = 24
dom.KernelType	= 0
dom.Scheme	    = 0
#     	dom.XSPH	= 0.5 //Very important

#double dx,h,rho,K,G,Cs,Fy
#double H,L,n

H	= 0.01
L	= 0.03
n	= 40.0

rho	= 1000.0
K	= 3.25e6
G	= 7.15e5
Fy	= 4000.0
dx	= H / n
h	= dx*1.3 #Very important
Cs	= sqrt(K/rho)

#double timestep
timestep = (0.2*h/(Cs))

# cout<<"t  = "<<timestep<<endl
# cout<<"Cs = "<<Cs<<endl
# cout<<"K  = "<<K<<endl
# cout<<"G  = "<<G<<endl
# cout<<"Fy = "<<Fy<<endl
#dom.GeneralAfter = & UserAcc
dom.DomMax[0] = L
dom.DomMin[0] = -L

dom.AddBoxLength(1 ,Vec3_t ( -L/2.0-L/20.0 , -H/2.0 , 0.0 ), L + L/10.0 + dx/10.0 , H + dx/10.0 ,  0 , dx/2.0 ,rho, h, 1 , 0 , False, False )

print ("Particle Number")
print (len(dom.Particles))
#double x

#for (size_t a=0 a<dom.Particles.Size() a++)
for a in range (len (dom.Particles)):
    dom.Particles[a].G		    = G
    dom.Particles[a].PresEq	    = 0
    dom.Particles[a].Cs		    = Cs
    dom.Particles[a].Shepard	= False
    dom.Particles[a].Material	= 2
    dom.Particles[a].Fail		= 1
    dom.Particles[a].Sigmay	    = Fy
    dom.Particles[a].Alpha		= 1.0
    dom.Particles[a].TI		    = 0.3
    dom.Particles[a].TIInitDist	= dx
    x = dom.Particles[a].x(0)
    if (x < -L/2.0):
        dom.Particles[a].ID=2
    if (x>L/2.0):
        dom.Particles[a].ID=3


#    	dom.WriteXDMF("maz")
#dom.Solve(/*tf*/1000.0,/*dt*/timestep,/*dtOut*/0.001,"test06",999)
dom.Solve(1000.0,timestep,0.001,"test06",999)
#return 0