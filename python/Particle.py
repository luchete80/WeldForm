from Vec3_t import *
from Functions import *
from numpy import *
import numpy.matlib

class Particle:
#inline Particle::Particle(int Tag, Vec3_t const & x0, Vec3_t const & v0, Mass0, Density0, h0,bool Fixed)

    def __init__(self, Tag, x0, v0, Mass0,  Density0,  h0, Fixed):
        self.ct = 0
        self.a = 0.0
        self.x = Vec3_t(x0)
        self.n = 0.0
        self.n0 = 0.0
        self.k = 0.0
        self.k2 = 0.0

        self.Cs		= 0.0
        self.P0		= 0.0
        self.PresEq	= 0 
        self.Alpha	= 0.0
        self.Beta	= 0.0

        self.Mu	=0.	#/< Dynamic viscosity coefficient of the fluid particle
        self.MuRef	=0.	#/< Reference Dynamic viscosity coefficient
        self.T0		=0.#/< Yield stress for Bingham fluids
        m		=0.#/< Normalization value for Bingham fluids
		# size_t	VisM		#/< Non-Newtonian viscosity method
		# bool		LES		#/< Large eddy simulation using sub-particle scale
		# double	CSmag		#/< Coefficient of Smagorinsky-Lilly model
        
        self.va = 0.0
        self.vb = 0.0
        self.NSv = 0.0
        self.FSINSv = 0.0
        self.v = v0
        self.VXSPH = 0.0
        self.TI		= 0.0
        self.TIn		= 4.0
        self.TIInitDist  = 0.0

        self.Densitya = 0.0
        self.Densityb = 0.0
        self.Density = Density0
        self.RefDensity = Density0

        self.Mass = Mass0
        self.FPMassC = 1.0
        #IsFree = !Fixed
        self.h = h0
        self.Pressure=0.0
        self.FSIPressure=0.0
        self.ID = Tag
        self.CC=zeros(3)
        self.CC[0]= self.CC[1] = self.CC[2] = 0
        self.LL=0
        self.ZWab = 0.0
        self.SumDen = 0.0
        self.dDensity=0.0
        self.ShearRate = 0.0
        self.MuRef = Mu = 0.0
        self.VisM = 0
        self.T0 = 0.0
        self.m = 300.0
        self.SumKernel = 0.0
        self.FSISumKernel = 0.0
        self.G = 0.0
        self.K = 0.0
        self.Material = 0
        self.Fail = 0
        self.c = 0.0
        self.phi = 0.0
        self.psi = 0.0
        self.d =0.0
        self.Sigmay = 0.0
        self.NoSlip = False
        self.Shepard = False
        self.InOut = 0
        self.FirstStep = True
        self.V = self.Mass/self.RefDensity
        self.RhoF = 0.0
        self.IsSat = False
        self.SatCheck = False
        self.ShepardStep = 40
        self.ShepardCounter = 0
        self.S = 0.0
        self.VarPorosity = False
        self.SeepageType = 0
        self.S = 0
        self.LES = False
        self.SBar = 0.0
        self.CSmag = 0.17
        
        ShearStress =matrix(numpy.matlib.zeros((3, 3)))
        StrainRate  =matrix(numpy.matlib.zeros((3, 3)))


        # set_to_zero(Strainb)
        # set_to_zero(Strain)
        # set_to_zero(Sigmab)
        # set_to_zero(Sigma)
        # #    set_to_zero(FSISigma)
        # set_to_zero(Sigmaa)
        # set_to_zero(ShearStress)
        # set_to_zero(ShearStressb)
        # set_to_zero(TIR)
        # set_to_zero(StrainRate)
        # set_to_zero(RotationRate)
        #omp_init_lock(&my_lock)
        
    #inline void Particle::Move(dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin, size_t Scheme, Mat3_t I)
    def Move(dt, Domainsize, domainmax, domainmin, Scheme, I):
        ds=Vec3_t(Domainsize)
        if (Scheme == 0):
            Move_MVerlet(I, dt)
        else:
            Move_Leapfrog(I, dt)


        #Periodic BC particle position update
        # if (Domainsize(0)>0.0):
            # x(0) = x(0) - Domainsize(0) if (x(0)>(domainmax(0))) else x(0)
            # (x(0)<(domainmin(0))) ? x(0) += Domainsize(0) : x(0)
        # if (Domainsize(1)>0.0):
            # (x(1)>(domainmax(1))) ? x(1) -= Domainsize(1) : x(1)
            # (x(1)<(domainmin(1))) ? x(1) += Domainsize(1) : x(1)
        # if (Domainsize(2)>0.0):
            # (x(2)>(domainmax(2))) ? x(2) -= Domainsize(2) : x(2)
            # (x(2)<(domainmin(2))) ? x(2) += Domainsize(2) : x(2)
        


    #inline void Particle::Move_MVerlet (Mat3_t I, dt)
    def Move_MVerlet (I, dt):
        if (FirstStep):
            ct = 30
            FirstStep = False

        x += dt*(v+VXSPH) + 0.5*dt*dt*a

        if (ct == 30):
            #if (Shepard==True && ShepardCounter == ShepardStep):
            if (Shepard==True and ShepardCounter == ShepardStep):
            
                if (ZWab > 0.6):                
                    Densityb	= SumDen/ZWab
    #				Densityb	= Density
                    Density		= SumDen/ZWab              
                else:                
                    Densityb	= Density
                    Density		+=dt*dDensity
            
            else:
            
                Densityb		= Density
                Density			+=dt*dDensity

            vb	= v
            v	+=dt*a
        
        else:
        
            if (Shepard==True and ShepardCounter == ShepardStep):
            
                if (ZWab>0.6):
                
                    Densityb	= SumDen/ZWab
    #				Densityb	= Density
                    Density		= SumDen/ZWab
                
                else:
                
                    dens	= Density
                    Density		= Densityb + 2.0*dt*dDensity
                    Densityb	= dens
                
            
            else:
            
                dens	= Density
                Density		= Densityb + 2.0*dt*dDensity
                Densityb	= dens
            

            #Vec3_t temp
            temp	= v
            v		= vb + 2*dt*a
            vb		= temp
        

        if Material==1:
            Mat1(dt)
        else:
            if Material == 2:
                Mat2MVerlet(dt)
            else:
                Mat3MVerlet(I,dt)
       # default:
            # std::cout << "Material Type No is out of range. Please correct it and run again" << std::endl
            # std::cout << "1 => Fluid" << std::endl
            # std::cout << "2 => Solid" << std::endl
            # std::cout << "3 => Soil" << std::endl
            # abort()
            # break
        
        if (ct == 30):
            ct = 0 
        else: 
            ct+=1
        if (ShepardCounter == ShepardStep): 
            ShepardCounter = 0 
        else: 
            ShepardCounter+=1
    
    #inline void Particle::Mat2MVerlet(dt)
    def Mat2MVerlet(dt):
    
        Pressure = EOS(PresEq, Cs, P0,Density, RefDensity)

        # Jaumann rate terms
        #Mat3_t RotationRateT, Stress,SRT,RS
        SRT=ShearStress*RotationRate.transpose()
        #Trans(RotationRate,RotationRateT)
        #Mult(ShearStress,RotationRateT,SRT)
        Mult(RotationRate,ShearStress,RS)

        # Elastic prediction step (ShearStress_e n+1)
        Stress			= ShearStress
        if (ct == 30):
            ShearStress	= dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*identity(3))+SRT+RS) + ShearStress
        else:
            ShearStress	= 2.0*dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*identity(3))+SRT+RS) + ShearStressb
        ShearStressb	= Stress

        if (Fail == 1):
            J2	= 0.5*(ShearStress(0,0)*ShearStress(0,0) + 2.0*ShearStress(0,1)*ShearStress(1,0) +
                            2.0*ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,1)*ShearStress(1,1) +
                            2.0*ShearStress(1,2)*ShearStress(2,1) + ShearStress(2,2)*ShearStress(2,2))
            #Scale back
            ShearStress	= min((Sigmay/sqrt(3.0*J2)),1.0)*ShearStress
        

        Sigma			= -Pressure * identity(3) + ShearStress

        Stress	= Strain
        if (ct == 30):
            Strain	= dt*StrainRate + Strain
        else:
            Strain	= 2.0*dt*StrainRate + Strainb
        Strainb	= Stress


        if (Fail > 1):      
            #std::cout<<"Undefined failure criteria for solids"<<std::endl
            print ("Undefined failure criteria for solids")
            #abort()
        
    


    #inline void Particle::translate(dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin)
    def translate(dt, Domainsize, domainmax, domainmin):
        x = x + dt*v + 0.5*dt*dt*a

        #Evolve velocity
        #Vec3_t temp
        temp = v
        v = vb + 2*dt*a
        vb = temp

        #Periodic BC particle position update
        #TODO: UNCOMMENT
        if (Domainsize(0)>0.0):
        
            if (x(0)>(domainmax(0))): 
                x[0]= x[0] - Domainsize[0] 
            else: x(0)
            if (x(0)<(domainmin(0))): 
                x[0] += Domainsize[0] 
            else: x(0)
        
        # if (Domainsize(1)>0.0):
        
            # (x(1)>(domainmax(1))) ? x(1) -= Domainsize(1) : x(1)
            # (x(1)<(domainmin(1))) ? x(1) += Domainsize(1) : x(1)
        
        # if (Domainsize(2)>0.0):
        
            # (x(2)>(domainmax(2))) ? x(2) -= Domainsize(2) : x(2)
            # (x(2)<(domainmin(2))) ? x(2) += Domainsize(2) : x(2)
        
    