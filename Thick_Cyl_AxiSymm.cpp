#include "Domain.h"
#include "InteractionAlt.cpp"
#include "SolverKickDrift.cpp"
#include "SolverFraser.cpp"

#define TAU		1.00
#define VMAX	0.05



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
		if (domi.Particles[i]->ID == 2)  //FIXED
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,-vcompress/2.);
			domi.Particles[i]->va		= Vec3_t(0.0,0.0,-vcompress/2.);
			domi.Particles[i]->vb		= Vec3_t(0.0,0.0,-vcompress/2.);

		}
    
	}
}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
		 SPH::Domain	dom;

		dom.Dimension	= 2;
    dom.dom_bid_type = AxiSymmetric ;
		dom.Nproc	= 12;
		dom.Kernel_Set(Qubic_Spline);
		dom.Scheme	= 1;	//Mod Verlet
		//dom.XSPH	= 0.1; //Very important

			double dx,h,rho,K,G,Cs,Fy;
		double R,L,n;

		////////Stable Axi Symmetric SPH formulation with no axis singularity
    ////////Jian Wang, Dave Chan, Jiahe Zhang
    float r1, r2;
    
    r1 = 0.2;
    r2 = 1.2;
    
		L	= 1.0;
		n	= 30.0;		//in length, radius is same distance
		
		rho	= 2700.0;
		K	= 6.7549e10;
		G	= 2.5902e10;
		Fy	= 300.e6;
    	//dx	= L / (n-1);
		//dx = L/(n-1);
		dx = 0.10;
    h	= dx*1.2; //Very important
        Cs	= sqrt(K/rho);

        double timestep;
        timestep = (0.7*h/(Cs));
		
		//timestep = 2.5e-6;

        cout<<"t  = "<<timestep<<endl;
        cout<<"Cs = "<<Cs<<endl;
        cout<<"K  = "<<K<<endl;
        cout<<"G  = "<<G<<endl;
        cout<<"Fy = "<<Fy<<endl;
    	dom.GeneralAfter = & UserAcc;
        dom.DomMax(0) = r2;
        dom.DomMin(0) = r1;
        dom.DomMax(1) = L;
        dom.DomMin(1) = 0.;
        
    dom.AddBoxLength(0, Vec3_t(r1,0.,0.), r2-r1, L, 0., 
									dx/2., rho, h,  1 , 0 , false, false );
                  
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

    		if ( x > r2 - 2.0*dx ){
    			dom.Particles[a]->ID=2;	         
        }

      }
		dom.WriteXDMF("maz");
		dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
		dom.BC.InOutFlow = 0;
    // timestep = (1.0*h/(Cs+VMAX)); 
    // dom.CFL = 1.0;
    // //timestep = 2.5e-6;
    // dom.auto_ts = false;
    
    //dom.Solve_orig_Ext(/*tf*/0.00205,/*dt*/timestep,/*dtOut*/0.001,"test06",999);
		//dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/0.0001,"test06",999);
    dom.SolveDiffUpdateFraser(/*tf*/0.105,/*dt*/timestep,/*dtOut*/1.e-2,"test06",1000);
    //dom.SolveDiffUpdateFraser(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/timestep,"test06",1000);
    
		return 0;
}
MECHSYS_CATCH
