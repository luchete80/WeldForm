#include "Domain.h"

#define TAU		0.5
#define UMAX	0.0024
 
double part_force,part_acc;
double accmax;

void UserAcc(SPH::Domain & domi)
{
	
	double acc;
	if (domi.getTime() < TAU ) 
		acc = accmax*domi.getTime()/TAU;
	else
		acc = accmax;
	
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)

	#ifdef __GNUC__
	for (size_t i=0; i<domi.Particles.Size(); i++)
	#else
	for (int i=0; i<domi.Particles.Size(); i++)
	#endif
	
	
	
	{
		if (domi.Particles[i]->ID == 3)
		{
			double Force = domi.Particles[i]->Mass;
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->vb		= Vec3_t(0.0,0.0,0.0);
//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
		if (domi.Particles[i]->ID == 2)
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->vb		= Vec3_t(0.0,0.0,0.0);
//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
		if (domi.Particles[i]->ID == 4)
		{
			domi.Particles[i]->a		= Vec3_t(0.0,-acc,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->vb		= Vec3_t(0.0,0.0,0.0);
//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
	}
}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
       SPH::Domain	dom;

        dom.Dimension	= 3;
        dom.Nproc	= 4;
    	dom.Kernel_Set(Quintic);
    	dom.Scheme	= 0;	//Mod Verlet
//     	dom.XSPH	= 0.5; //Very important
		
        double dx,h,rho,K,G,Cs,Fy;
    	double H,Lx,Ly,Lz;
		
		Lx = 10.0;
		Ly = 1.;
		Lz = 1.;	
		
		double E  = 70.e9;
		double nu = 0.33;
		
    	rho	= 2700.0;
		K= E / ( 3.*(1.-2*nu) );
		G= E / (2.* (1.+nu));
		Fy	= 1000.e10;

		dx = 0.111;
    	h	= dx*1.1; //Very important
        Cs	= sqrt(K/rho);

        double timestep;
        timestep = (1.*h/(Cs));
		
		//timestep = 2.5e-6;
		//timestep = 5.e-7;

        cout<<"t  = "<<timestep<<endl;
        cout<<"Cs = "<<Cs<<endl;
        cout<<"K  = "<<K<<endl;
        cout<<"G  = "<<G<<endl;
        cout<<"Fy = "<<Fy<<endl;
    	dom.GeneralAfter = & UserAcc;
        dom.DomMax(0) = Lx;
        dom.DomMin(0) = -Lx;

		// inline void Domain::AddBoxLength(int tag, Vec3_t const & V, 
									//double Lx, double Ly, double Lz, 
									// double r, double Density, double h, int type, int rotation, bool random, bool Fixed)
									

     	dom.AddBoxLength(1 ,Vec3_t ( -Lx/2.0 -Lx/40.0, -Ly/2.0 , -Lz/2.0 ), 
							Lx + Lx/20.0 /*+ dx*/, Ly /*+dx*/,  Lz /*+dx*/, 
							dx/2.0 ,rho, h, 1 , 0 , false, false );


		cout << "Particle count: "<<dom.Particles.Size()<<endl;
			int load_part_count = 0;
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
    		dom.Particles[a]->TI		= 0.3;
    		dom.Particles[a]->TIInitDist	= dx;
    		double x = dom.Particles[a]->x(0);
			double y = dom.Particles[a]->x(1);
			double z = dom.Particles[a]->x(2);
			//cout << "particle Mass"<<dom.Particles[a]->Mass<<endl;

    		if ( x <= -Lx/2. )
    			dom.Particles[a]->ID=2;
    		else if ( /*y >= (Ly/2. -3*dx ) && */x >= Lx/2.  )
    			dom.Particles[a]->ID=3;
			else if (y > Ly/2.0-dx){
				load_part_count++;
				dom.Particles[a]->ID=4;
			}

    	}
		cout << "Load Part Count: "<<load_part_count<<endl;		
		part_force = 500.e3/load_part_count;
		accmax = part_force/dom.Particles[0]->Mass;
		cout << "Particle Force" << part_force<<"Max Acc: "<<accmax<<endl;
		
		dom.WriteXDMF("maz");
//		dom.m_kernel = SPH::iKernel(dom.Dimension,h);	


    	dom.Solve(/*tf*/1.01,/*dt*/timestep,/*dtOut*/0.001,"test06",999);
        return 0;
}
MECHSYS_CATCH
