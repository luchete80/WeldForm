from numpy import *
from Vec3_t import *
from Particle import *
from Boundary import *

#https:#stackoverflow.com/questions/32558805/ceil-and-floor-equivalent-in-python-3-without-math-module
def ceil(n):
    res = int(n)
    return res if res == n or n < 0 else res+1

def floor(n):
    res = int(n)
    return res if res == n or n >= 0 else res-1


# def ceil(n):
    # return int(-1 * n # 1 * -1)

# def floor(n):
    # return int(n # 1)   
    
class Domain:
    
    def __init__(self):
        self.OutputName=[0 for i in range(3)] 
        self.OutputName[0] = "Property1"
        self.OutputName[1] = "Property2"
        self.OutputName[2] = "Property3"
        
        self.Time    = 0.0

        self.Dimension = 2


        self.Gravity	= Vec3_t(0.0,0.0,0.0)

        self.MuMax=0.		#/< Max Dynamic viscosity for calculating the timestep
        self.CsMax=0.		#/< Max speed of sound for calculating the timestep
    
        self.Cellfac = 2.0

        self.KernelType	= 0
        self.SWIType	= 0
        #self.FSI        = false
        self.VisEq	    = 0
        self.Scheme	= 0
        self.GradientType = 0  

        XSPH	= 0.0
        InitialDist = 0.0

        AvgVelocity = 0.0


        #omp_init_lock (&dom_lock)
        Nproc	= 1
        
        self.I=matrix(numpy.matlib.zeros((3, 3)))

        self.deltat	= 0.0
        self.deltatint	= 0.0
        self.deltatmin	= 0.0
        self.sqrt_h_a = 0.0025

        self.TRPR = Vec3_t(0.,0.,0.)
        self.BLPF = Vec3_t(0.,0.,0.)

        self.CellSize=Vec3_t(0.,0.,0.);      	#/< Calculated cell size according to (cell size >= 2h)
        self.CellNo=[0,0,0];      #/< No. of cells for linked list
        self.hmax=0.;		#/< Max of h for the cell size  determination
        self.DomSize	= Vec3_t(0.0,0.0,0.0)
        self.rhomax=0.;        
        
        self.DomMax = Vec3_t(-100000000000,-100000000000,-100000000000)
        self.DomMin = Vec3_t(100000000000,100000000000,100000000000)
        
        x=Vec3_t(0,0,0)
        print(x)
        self.Particles=[]
        self.AddSingleParticle(1,x,1,1,1)
        
        self.BC=Boundary()
        
    def AddSingleParticle(tag, x, Mass, Density, h, Fixed):
    #int tag, Vec3_t const & x, double Mass, double Density, double h, bool Fixed)
        print(x)
        print("Hello")
    
    def myfunc(self):
        print("Hello my name is " + self.name)
        

    #inline void Domain::AddSingleParticle(int tag, Vec3_t const & x, double Mass, double Density, double h, bool Fixed)
    
        #Particles.Push(new Particle(tag,x,Vec3_t(0,0,0),Mass,Density,h,Fixed))
            
    #inline void Domain::AddBoxLength(int tag, Vec3_t const & V, double Lx, double Ly, double Lz, double r, double Density, double h, int type, int rotation, bool random, bool Fixed)
    def AddBoxLength(self,tag, V, Lx, Ly, Lz, r, Density, h, type, rotation, random, Fixed):
        if (not(type==0 or type==1)):
        
            # print( "Packing Type is out of range. Please correct it and run again" )
            # print( "0 => Hexagonal Close Packing" )
            # print( "1 => Cubic Packing" )
            return
        

        if (not(rotation==0 or rotation==90)):
        
            print( "Packing Rotation Angle is out of range. Please correct it and run again" )
            print( "0 => " )
            print( "0 0 0 0" )
            print( " 0 0 0 0" )
            print( "0 0 0 0" )
            print( " 0 0 0 0" )
            print( "90 => Cubic Close Packing" )
            print( "  0   0" )
            print( "0 0 0 0" )
            print( "0 0 0 0" )
            print( "0   0  " )
            abort()
        

    #	Util::Stopwatch stopwatch
        print( "\n--------------Generating particles by AddBoxLength with defined length of particles-----------" )

        #size_t PrePS = Particles.Size()
        #PrePS=len(Particles)   #TODO: CHECK DIMENSION
        PrePS=0
        
        #double x,y,xp,yp
        #size_t i,j

        #double qin = 0.03
        qin = 0.03
        #srand(100)
        #random(100)

        if (self.Dimension==3):        
            if (type==0):
            
                #Hexagonal close packing
                #double z,zp
                #size_t k=0
                k=0
                zp = V(2)   #TODO: CHECK WHAT TO DO!

                while (zp <= (V(2)+Lz-r)):
                
                    j = 0
                    yp = V(1)
                    while (yp <= (V(1)+Ly-r)):
                    
                        i = 0
                        xp = V(0)
                        while (xp <= (V(0)+Lx-r)):
                        
                            
                            if ((k%2!=0) and (j%2!=0)) :
                                x = V(0) + (2*i+(j%2)+(k%2)-1)*r 
                            else: 
                                x = V(0) + (2*i+(j%2)+(k%2)+1)*r
                            y = V(1) + (sqrt(3.0)*(j+(1.0/3.0)*(k%2))+1)*r
                            z = V(2) + ((2*sqrt(6.0)/3)*k+1)*r
                            #TODO: UNCOMMENT BELOW!
                            # if (random):
                                # Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed))
                            # else:    
                                # Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),0.0,Density,h,Fixed))
                            i+=1
                            if ((k%2!=0) and (j%2!=0)):
                                xp = V(0) + (2*i+(j%2)+(k%2)-1)*r 
                            else: 
                                xp = V(0) + (2*i+(j%2)+(k%2)+1)*r
                        
                        j+=1
                        yp = V(1) + (sqrt(3.0)*(j+(1.0/3.0)*(k%2))+1)*r
                    
                    k+=1
                    zp = V(2) + ((2*sqrt(6.0)/3)*k+1)*r
                
            
            else:
            
                #Cubic packing
                #double z,zp
                #size_t k=0
                k=0
                zp = V(2)

                while (zp <= (V(2)+Lz-r)):
                
                    j = 0
                    yp = V(1)
                    while (yp <= (V(1)+Ly-r)):
                    
                        i = 0
                        xp = V(0)
                        while (xp <= (V(0)+Lx-r)):
                        
                            x = V(0) + (2.0*i+1)*r
                            y = V(1) + (2.0*j+1)*r
                            z = V(2) + (2.0*k+1)*r
                            if (random):
                                x1=x + qin*r*double(rand())/RAND_MAX
                                y1=y+ qin*r*(rand())/RAND_MAX
                                z1=z+ qin*r*double(rand())/RAND_MAX
                                part=Particle(tag,Vec3_t(x1,y1,z1),Vec3_t(0,0,0),0.0,Density,h,Fixed)
                                Particles.append(part)
                                #Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed))
                            else:    
                                part=Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),0.0,Density,h,Fixed)
                                Particles.append(part)
                                #Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),0.0,Density,h,Fixed))
                            i+=1
                            xp = V(0) + (2*i+1)*r
                        
                        j+=1
                        yp = V(1) + (2.0*j+1)*r
                    
                    k+=1
                    zp = V(2) + (2.0*k+1)*r
                
            

            #Calculate particles' mass in 3D
            #Vec3_t temp, Max=V
            Max=V
            #for (size_t i=PrePS i<Particles.Size() i++)
            
            # for i in range(PrePS,Size(Particles)):
                # if (Particles[i].x(0) > Max(0)): Max(0) = Particles[i].x(0)
                # if (Particles[i].x(1) > Max(1)): Max(1) = Particles[i].x(1)
                # if (Particles[i].x(2) > Max(2)): Max(2) = Particles[i].x(2)
            
            Max +=r
            temp = Max-V
            #double Mass = temp(0)*temp(1)*temp(2)*Density/(Particles.Size()-PrePS)
            Mass = temp(0)*temp(1)*temp(2)*Density/(Size(Particles)-PrePS)

            #pragma omp parallel for num_threads(Nproc)
            #for (size_t i=PrePS i<Particles.Size() i++)
            for i in range(PrePS,Size(Particles)):
                Particles[i].Mass = Mass
            
        

        if (self.Dimension==2):
        
            if (type==0):
            
                #Hexagonal close packing
                if (rotation==0):
                
                    j = 0
                    yp = V(1)

                    while (yp <= (V(1)+Ly-r)):
                    
                        i = 0
                        xp = V(0)
                        while (xp <= (V(0)+Lx-r)):
                        
                            x = V(0) + (2*i+(j%2)+1)*r
                            y = V(1) + (sqrt(3.0)*j+1)*r
                            if (random):
                                x1=x + qin*r*double(rand())/RAND_MAX
                                y1=y+ qin*r*(rand())/RAND_MAX
                                part=Particle(tag,Vec3_t(x1,y1,0.),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed)
                                #Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed))
                                Particles.append(part)
                            else: 
                                part=Particle(tag,Vec3_t(x,y,0.),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed)
                                Particles.append(part)
                                #Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed))
                            i+=1
                            xp = V(0) + (2*i+(j%2)+1)*r
                        
                        j+=1
                        yp = V(1) + (sqrt(3.0)*j+1)*r
                    
                
                else:
                
                    i = 0
                    xp = V(0)

                    while (xp <= (V(0)+Lx-r)):
                    
                        j = 0
                        yp = V(1)
                        while (yp <= (V(1)+Ly-r)):
                        
                            x = V(0) + (sqrt(3.0)*i+1)*r
                            y = V(1) + (2*j+(i%2)+1)*r
                            if (random): 
                                x1=x + qin*r*double(rand())/RAND_MAX
                                y1=y+ qin*r*(rand())/RAND_MAX
                                part=Particle(tag,Vec3_t(x1,y1,0.),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed)
                                Particles.append(part)
                                #Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed))
                            else:    
                                #Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed))
                                part=Particle(tag,Vec3_t(x,y,0.),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed)
                                Particles.append(part)
                            j+=1
                            yp = V(1) + (2*j+(i%2)+1)*r
                        
                        i+=1
                        xp = V(0) + (sqrt(3.0)*i+1)*r
                    
                
            
            else:
            
                #Cubic packing
                j = 0
                yp = V[1]

                while (yp <= (V[1]+Ly-r)):
                
                    i = 0
                    xp = V(0)
                    while (xp <= (V(0)+Lx-r)):
                    
                        x = V(0) + (2*i+1)*r
                        y = V(1) + (2*j+1)*r
                        if (random) :
                            x1=x + qin*r*double(rand())/RAND_MAX
                            y1=y+ qin*r*(rand())/RAND_MAX
                            part=Particle(tag,Vec3_t(x1,y1,0.),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed)
                            self.Particles.append(part)
                            #Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed))
                        else:    
                            part=Particle(tag,Vec3_t(x,y,0.),Vec3_t(0,0,0),2.0*r*2.0*r*Density,Density,h,Fixed)
                            self.Particles.append(part)
                            #Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),2.0*r*2.0*r*Density,Density,h,Fixed))
                        i+=1
                        xp = V(0) + (2*i+1)*r
                    
                    j+=1
                    yp = V(1) + (2*j+1)*r
                

            
        

        R = r
        
        
    #inline void Domain::CellInitiate ()
    def CellInitiate(self):
            
        hola=float(self.TRPR*self.TRPR)**0.5
        print(self.TRPR)
        self.TRPR.norm()
        
        if (not(float(self.TRPR.norm())>0.0) and not(float(self.BLPF.norm())>0.0)):
        
            # Calculate Domain Size
            self.BLPF = self.Particles[0].x
            self.TRPR = self.Particles[0].x
            self.hmax = self.Particles[0].h
            self.rhomax = self.Particles[0].Density

            #for (size_t i=0 i<self.Particles.Size() i++)
            for i in range(len(self.Particles)):
            
                if (self.Particles[i].x(0) > self.TRPR(0)): self.TRPR[0] = self.Particles[i].x(0)
                if (self.Particles[i].x(1) > self.TRPR(1)): self.TRPR[1] = self.Particles[i].x(1)
                if (self.Particles[i].x(2) > self.TRPR(2)): self.TRPR[2] = self.Particles[i].x(2)

                if (self.Particles[i].x(0) < self.BLPF(0)): self.BLPF[0] = self.Particles[i].x(0)
                if (self.Particles[i].x(1) < self.BLPF(1)): self.BLPF[1] = self.Particles[i].x(1)
                if (self.Particles[i].x(2) < self.BLPF(2)): self.BLPF[2] = self.Particles[i].x(2)

                if (self.Particles[i].h > self.hmax): self.hmax=self.Particles[i].h
                if (self.Particles[i].Density > self.rhomax): self.rhomax=self.Particles[i].Density
                if (self.Particles[i].Mu > self.MuMax): self.MuMax=self.Particles[i].Mu
                if (self.Particles[i].Cs > self.CsMax): self.CsMax=self.Particles[i].Cs
            
        

        # Override the calculated domain size
        if (self.DomMax(0)>self.TRPR(0)): self.TRPR[0] = self.DomMax(0)
        if (self.DomMax(1)>self.TRPR(1)): self.TRPR[1] = self.DomMax(1)
        if (self.DomMax(2)>self.TRPR(2)): self.TRPR[2] = self.DomMax(2)
        if (self.DomMin(0)<self.BLPF(0)): self.BLPF[0] = self.DomMin(0)
        if (self.DomMin(1)<self.BLPF(1)): self.BLPF[1] = self.DomMin(1)
        if (self.DomMin(2)<self.BLPF(2)): self.BLPF[2] = self.DomMin(2)


        #Because of Hexagonal close packing in x direction domain is modified
        #TODO:UNCOMMENT BELOW!
        # if (!self.BC.Periodic[0]): self.TRPR(0) += self.hmax/2	self.BLPF(0) -= self.hmax/2 else: self.TRPR(0) += R self.BLPF(0) -= R
        # if (!self.BC.Periodic[1]): self.TRPR(1) += self.hmax/2	self.BLPF(1) -= self.hmax/2 else: self.TRPR(1) += R self.BLPF(1) -= R
        # if (!self.BC.Periodic[2]): self.TRPR(2) += self.hmax/2	self.BLPF(2) -= self.hmax/2 else: self.TRPR(2) += R self.BLPF(2) -= R

        # Calculate Cells Properties
        if (self.Dimension == 2):
            if (double (ceil(((self.TRPR(0)-self.BLPF(0))/(self.Cellfac*self.hmax)))-((self.TRPR(0)-self.BLPF(0))/(self.Cellfac*self.hmax)))<(self.hmax/10.0)):
                self.CellNo[0] = int(ceil((self.TRPR(0)-self.BLPF(0))/(self.Cellfac*self.hmax)))
            else:
                self.CellNo[0] = int(floor((self.TRPR(0)-self.BLPF(0))/(self.Cellfac*self.hmax)))

            if (float (ceil(((self.TRPR(1)-self.BLPF(1))/(self.Cellfac*self.hmax)))-((self.TRPR(1)-self.BLPF(1))/(self.Cellfac*self.hmax)))<(self.hmax/10.0)):
                self.CellNo[1] = int(ceil((self.TRPR(1)-self.BLPF(1))/(self.Cellfac*self.hmax)))
            else:
                self.CellNo[1] = int(floor((self.TRPR(1)-self.BLPF(1))/(self.Cellfac*self.hmax)))

            self.CellNo[2] = 1
            
            print("self.CellNo[0]",self.CellNo[0])#LUCIANO
            CellSize  = Vec3_t ((self.TRPR(0)-self.BLPF(0))/self.CellNo[0],(self.TRPR(1)-self.BLPF(1))/self.CellNo[1],0.0)

        elif (self.Dimension == 3):
            if (double (ceil(((self.TRPR(0)-self.BLPF(0))/(self.Cellfac*self.hmax)))-((self.TRPR(0)-self.BLPF(0))/(self.Cellfac*self.hmax)))<(self.hmax/10.0)):
                self.CellNo[0] = int(ceil((self.TRPR(0)-self.BLPF(0))/(self.Cellfac*self.hmax)))
            else:
                self.CellNo[0] = int(floor((self.TRPR(0)-self.BLPF(0))/(self.Cellfac*self.hmax)))

            if (double (ceil(((self.TRPR(1)-self.BLPF(1))/(self.Cellfac*self.hmax)))-((self.TRPR(1)-self.BLPF(1))/(self.Cellfac*self.hmax)))<(self.hmax/10.0)):
                self.CellNo[1] = int(ceil((self.TRPR(1)-self.BLPF(1))/(self.Cellfac*self.hmax)))
            else:
                self.CellNo[1] = int(floor((self.TRPR(1)-self.BLPF(1))/(self.Cellfac*self.hmax)))

            if (double (ceil(((self.TRPR(2)-self.BLPF(2))/(self.Cellfac*self.hmax)))-((self.TRPR(2)-self.BLPF(2))/(self.Cellfac*self.hmax)))<(self.hmax/10.0)):
                self.CellNo[2] = int(ceil((self.TRPR(2)-self.BLPF(2))/(self.Cellfac*self.hmax)))
            else:
                self.CellNo[2] = int(floor((self.TRPR(2)-self.BLPF(2))/(self.Cellfac*self.hmax)))

            CellSize  = Vec3_t ((self.TRPR(0)-self.BLPF(0))/self.CellNo[0],(self.TRPR(1)-self.BLPF(1))/self.CellNo[1],(self.TRPR(2)-self.BLPF(2))/self.CellNo[2])

        else:
            print("Please correct the dimension (2=>2D or 3=>3D) and run again")
            #abort()

        

        # Periodic BC modifications
        if (BC.Periodic[0]): self.CellNo[0] += 2
        if (BC.Periodic[1]): self.CellNo[1] += 2
        if (BC.Periodic[2]): self.CellNo[2] += 2

        if (BC.Periodic[0]): DomSize[0] = (self.TRPR(0)-self.BLPF(0))
        if (BC.Periodic[1]): DomSize[1] = (self.TRPR(1)-self.BLPF(1))
        if (BC.Periodic[2]): DomSize[2] = (self.TRPR(2)-self.BLPF(2))

        # Initiate Head of Chain array for Linked-List
        #TODO: UNCOMMENT BELOW!
        # HOC = new int**[(int) self.CellNo[0]]
        # for(int i =0 i<self.CellNo[0] i++)
           # self.HOC[i] = new int*[self.CellNo[1]]
           # for(int j =0 j<self.CellNo[1] j++)
               # self.HOC[i][j] = new int[self.CellNo[2]]
               # for(int k = 0 k<self.CellNo[2]k++)
                  # self.HOC[i][j][k] = -1
               
           
        
        # Initiate Pairs array for neibour searching
        #TODO: UNCOMMENT BELOW!
        # for(size_t i=0  i<Nproc  i++)
        
        # self.SMPairs.Push(Initial)
        # self.NSMPairs.Push(Initial)
        # self.FSMPairs.Push(Initial)
        
    
    
    #inline void Domain::InitialChecks()
    def InitialChecks(self):
        #initializing identity matrix
        if (self.Dimension == 2): self.I[2,2] = 0

        if (self.Dimension<=1 or self.Dimension>3):
        
            print( "Please correct the dimension (2=>2D or 3=>3D) and run again" )
            #abort()
        

        if (self.BC.InOutFlow>0 and self.BC.Periodic[0]):
            print("FATAL: Periodic BC in the X direction cannot be used with In/Out-Flow BC simultaneously")

        #pragma omp parallel for schedule (static) num_threads(Nproc)
        #for (size_t i=0 i<self.Particles.Size() i++)
        for i in range(len(self.Particles)):
        
            #Initializing pressure of solid and fluid self.Particles
            if (self.Particles[i].Material < 3):
                self.Particles[i].Pressure = EOS(self.Particles[i].PresEq, self.Particles[i].Cs, self.Particles[i].P0,self.Particles[i].Density, self.Particles[i].RefDensity)

            # Initializing the permeability for soil self.Particles
            if (self.Particles[i].Material == 3):
            
                if(self.Particles[i].SeepageType == 1):
                    self.Particles[i].k = self.Particles[i].n0*self.Particles[i].n0*self.Particles[i].n0*self.Particles[i].d*self.Particles[i].d/(180.0*(1.0-self.Particles[i].n0)*(1.0-self.Particles[i].n0))

                elif self.Particles[i].SeepageType == 2:
                    self.Particles[i].k = self.Particles[i].n0*self.Particles[i].n0*self.Particles[i].n0*self.Particles[i].d*self.Particles[i].d/(150.0*(1.0-self.Particles[i].n0)*(1.0-self.Particles[i].n0))
                    self.Particles[i].k2= 1.75*(1.0-self.Particles[i].n0)/(self.Particles[i].n0*self.Particles[i].n0*self.Particles[i].n0*self.Particles[i].d)


                elif self.Particles[i].SeepageType == 3:
                    self.Particles[i].k = self.Particles[i].n0*self.Particles[i].n0*self.Particles[i].n0*self.Particles[i].d*self.Particles[i].d/(150.0*(1.0-self.Particles[i].n0)*(1.0-self.Particles[i].n0))
                    self.Particles[i].k2= 0.4/(self.Particles[i].n0*self.Particles[i].n0*self.Particles[i].d)

                else:
                    print( "Seepage Type No is out of range. Please correct it and run again" )
                    print( "0 => Darcy's Law" )
                    print( "1 => Darcy's Law & Kozeny–Carman Eq" )
                    print( "2 => The Forchheimer Eq & Ergun Coeffs" )
                    print( "3 => The Forchheimer Eq & Den Adel Coeffs" )
                    abort()
                        
                

                self.Particles[i].n = self.Particles[i].n0
            
        
    

    #inline void Domain::Solve (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx)
    def Solve (self, tf, dt, dtOut, TheFileKey, maxidx):
    
        print( "\n--------------Solving---------------------------------------------------------------")

        idx_out = 1
        tout = self.Time

        #Initializing adaptive time step variables
        self.deltat = self.deltatint = self.deltatmin	= dt
    
    
        #TODO: UNCOMMENT ALL BELOW!!!!!
        self.InitialChecks()
        self.CellInitiate()
        # ListGenerate()
        # PrintInput(TheFileKey)
        # TimestepCheck()
        # WholeVelocity()


        # #Initial model output
        # if (TheFileKey!=NULL)
        
            # String fn
            # fn.Printf    ("%s_Initial", TheFileKey)
            # WriteXDMF    (fn.CStr())
            # print( "\nInitial Condition has been generated\n")
        

        # while (Time<tf && idx_out<=maxidx):
        
            # StartAcceleration(Gravity)
            # if (BC.InOutFlow>0) InFlowBCFresh()
            # MainNeighbourSearch()
            # GeneralBefore(*this)
            # PrimaryComputeAcceleration()
            # LastComputeAcceleration()
            # GeneralAfter(*this)

            # # output
            # if (Time>=tout):
            
                # if (TheFileKey!=NULL):
                
                    # String fn
                    # fn.Printf    ("%s_%04d", TheFileKey, idx_out)
                    # WriteXDMF    (fn.CStr())
                    # print( "\nOutput No. " << idx_out << " at " << Time << " has been generated")
                    # print( "Current Time Step = " <<deltat<<std::endl
                
                # idx_out++
                # tout += dtOut
            

            # AdaptiveTimeStep()
            # Move(deltat)
            # Time += deltat
            # if (BC.InOutFlow>0) InFlowBCLeave() else CheckParticleLeave ()
            # CellReset()
            # ListGenerate()
        

        # print( "\n--------------Solving is finished---------------------------------------------------")

    