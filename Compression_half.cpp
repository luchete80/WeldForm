

#include "Domain.h"

#define TAU		0.005
#define VMAX	10.0



void UserAcc(SPH::Domain & domi)
{
	double vcompress;

	if (domi.getTime() < TAU ) 
		vcompress = VMAX/TAU * domi.getTime();
	else
		vcompress = VMAX;
	//cout << "time: "<< domi.getTime() << "V compress "<< vcompress <<endl;
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)

	#ifdef __GNUC__
	for (size_t i=0; i<domi.Particles.Size(); i++)
	#else
	for (int i=0; i<domi.Particles.Size(); i++)
	#endif
	
	{
    // TOP
		if (domi.Particles[i]->ID == 11)  //FIXED
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,-vcompress/2.);
			domi.Particles[i]->va		= Vec3_t(0.0,0.0,-vcompress/2.);
			domi.Particles[i]->vb		= Vec3_t(0.0,0.0,-vcompress/2.);

		}


    // if (domi.Particles[i]->ID == 3) { //xy
			// domi.Particles[i]->a[2]		= 0.0; 
      // domi.Particles[i]->v[2] = domi.Particles[i]->va[2] = domi.Particles[i]->vb[2]		= 0.;
		// } 


    
	}
}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
		 SPH::Domain	dom;

		dom.Dimension	= 3;
		dom.Nproc	= 4;
		dom.Kernel_Set(Qubic_Spline);
		dom.Scheme	= 1;	//Mod Verlet
		//dom.XSPH	= 0.1; //Very important

			double dx,h,rho,K,G,Cs,Fy;
		double R,L,n;

		R	= 0.15;
		L	= 0.56;
		n	= 30.0;		//in length, radius is same distance
		
		rho	= 2700.0;
		K	= 6.7549e10;
		G	= 2.5902e10;
		Fy	= 300.e6;
    	//dx	= L / (n-1);
		//dx = L/(n-1);
		dx = 0.015;
    h	= dx*1.2; //Very important
        Cs	= sqrt(K/rho);

        double timestep;
        timestep = (0.2*h/(Cs));
		
		//timestep = 2.5e-6;

        cout<<"t  = "<<timestep<<endl;
        cout<<"Cs = "<<Cs<<endl;
        cout<<"K  = "<<K<<endl;
        cout<<"G  = "<<G<<endl;
        cout<<"Fy = "<<Fy<<endl;
    	dom.GeneralAfter = & UserAcc;
        dom.DomMax(0) = L;
        dom.DomMin(0) = -L;


		// inline void Domain::AddCylinderLength(int tag, Vec3_t const & V, double Rxy, double Lz, 
									// double r, double Density, double h, bool Fixed) {
		bool symlength = true; //Also Symmetric on z axis
		bool Fixed = false;
    bool ghost = false;
	//Cylinder Slice of 90 degree and half length
	// THIS DOES NOT HAVE Z INITIAL POSITION SINCE IT IS ZERO OR -LZ ACCORGIND TO LAST ARGUMENT
	// void AddDoubleSymCylinderLength(int tag, double Rxy, double Lz, 
									// double r, double Density, double h, bool Fixed, bool symlength = false);
									
		//dom.AddDoubleSymCylinderLength(1, R, L/2. + L/18.,  dx/2., rho, h, Fixed, symlength); 
		dom.AddCylinderLength(1, Vec3_t(0.,0.,0.), R, L/2. + L/15.,  dx/2., rho, h, Fixed, ghost); 
    
    dom.gradKernelCorr = false;
        
		cout << "Particle count: "<<dom.Particles.Size()<<endl;

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		dom.Particles[a]->G		= G;
    		dom.Particles[a]->PresEq	= 0;
    		dom.Particles[a]->Cs		= Cs;
    		dom.Particles[a]->Shepard	= false;
    		dom.Particles[a]->Material	= 2;
    		dom.Particles[a]->Fail		= 1;
    		dom.Particles[a]->Sigmay	= Fy;
    		dom.Particles[a]->Alpha		= 1.0;
    		//dom.Particles[a]->Beta		= 1.0;
    		dom.Particles[a]->TI		= 0.3;
    		dom.Particles[a]->TIInitDist	= dx;
    		double x = dom.Particles[a]->x(0);
    		double y = dom.Particles[a]->x(1);
    		double z = dom.Particles[a]->x(2);
        
        
        //BOTTOM PLANE
        if ( z < dx  && z > -dx/2. )
    			dom.Particles[a]->ID=3;
    		
        // if ( x < dx  && x > -dx/2. && z < L/2. - dx)
    			// dom.Particles[a]->ID=1;
    		// if ( y < dx  && y > -dx/2. && z < L/2. - dx)
    			// dom.Particles[a]->ID=2; 
          
        //x,y, central symmetry
    		// if ( y < dx  && y > -dx/2. && x < dx  && x > -dx/2. && z < L/2. - dx)
    			// dom.Particles[a]->ID=4;  
        
    		// if ( y < dx  && y > -dx/2. && z < dx  && z > -dx/2. ) //yz -5
    			// dom.Particles[a]->ID=5;           
     		// if (  x < dx  && x > -dx/2. && z < dx  && z > -dx/2. ) //xz - 6
    			// dom.Particles[a]->ID=6;         
        // if ( y < dx  && y > -dx/2. && x < dx  && x > -dx/2. && z < dx  && z > -dx/2. ) //xz - 7
    			// dom.Particles[a]->ID=7;   


    		if ( z > L/2. - dx ){
    			dom.Particles[a]->ID=11;	         
        }
        //TOP
    		// if ( y < dx  && y > -dx/2. && z > L/2. - dx ) //yz -5
    			// dom.Particles[a]->ID=8;           
     		// if (  x < dx  && x > -dx/2. && z > L/2. - dx ) //xz - 6
    			// dom.Particles[a]->ID=9;         
        // if ( y < dx  && y > -dx/2. && x < dx  && x > -dx/2. && z > L/2. - dx ) //xyz - 7
    			// dom.Particles[a]->ID=10;   	

      }
		dom.WriteXDMF("maz");
		dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
		dom.BC.InOutFlow = 0;

    //dom.Solve_orig_Ext(/*tf*/0.00205,/*dt*/timestep,/*dtOut*/0.001,"test06",999);
		dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/0.0001,"test06",999);
    
		return 0;
}
MECHSYS_CATCH
