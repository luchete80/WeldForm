#include "Mesh.h"
#include "Domain.h"
#include <iostream>

#define TAU		0.005
#define VMAX	0.1

using namespace SPH;
using namespace std;

void UserAcc(SPH::Domain & domi) {
	double vcompress;
		vcompress = VMAX;

  int first = 3588;
  int last = 3863;
  double dS = 0.0008*0.0008;
  domi.m_scalar_prop = 0.;
  for (int i = first;i<=last;i++){
    domi.m_scalar_prop += domi.Particles[i]->Sigma (2,2) * dS;
  }
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	#ifdef __GNUC__
	for (size_t i=0; i<domi.Particles.Size(); i++)
	#else
	for (int i=0; i<domi.Particles.Size(); i++)
	#endif	
	{
		if (domi.Particles[i]->ID == 11)  //FIXED
		{
      domi.Particles[i]->a		= Vec3_t(0.,0.,0.);
      domi.Particles[i]->v		= Vec3_t(0.,0.,-vcompress);
			domi.Particles[i]->va		= Vec3_t(0.,0.,-vcompress);
			domi.Particles[i]->vb		= Vec3_t(0.,0.,-vcompress);

		}
  }
}


int main() try{
	//
	TriMesh mesh;

	cout << "Creating Mesh" << endl;


	 SPH::Domain	dom;

	dom.Dimension	= 3;
	dom.Nproc	= 4;
	dom.Kernel_Set(Qubic_Spline);
	dom.Scheme	= 1;	//Mod Verlet
	//dom.XSPH	= 0.1; //Very important

		double dx,h,rho,K,G,Cs,Fy;
	double R,L,n;

	R	= 0.008;
	L	= 0.02185;
	n	= 30.0;		//in length, radius is same distance

	rho	= 7850.0;
  double E = 200.e9;
  double nu = 0.3;
  K= E / ( 3.*(1.-2*nu) );
  G= E / (2.* (1.+nu));
	Fy	= 900.e6;
	//dx	= L / (n-1);
	//dx = L/(n-1);
	dx = 0.0008;  //Tenth of radius
	h	= dx*1.2; //Very important
	Cs	= sqrt(K/rho);

	double timestep;
	timestep = (0.05*h/(Cs+VMAX)); //CHANGED WITH VELOCITY

//timestep = 2.5e-6;

	cout<<"t  = "<<timestep<<endl;
	cout<<"Cs = "<<Cs<<endl;
	cout<<"K  = "<<K<<endl;
	cout<<"G  = "<<G<<endl;
	cout<<"Fy = "<<Fy<<endl;
	dom.GeneralAfter = & UserAcc;
	dom.DomMax(0) = L;
	dom.DomMin(0) = -L;

	bool ghost = true;								
	dom.AddCylinderLength(0, Vec3_t(0.,0.,0.), R, L/2.,  dx/2., rho, h, false, ghost); 
	cout << "Max z plane position: " <<dom.Particles[dom.Particles.Size()-1]->x(2)<<endl;

	
	dom.ts_nb_inc = 5;
	dom.gradKernelCorr = true;
  
  dom.Particles[3863]->print_history = true;
  
  cout << "Particle 3589 coords xyz "<<dom.Particles[3863]->x<<endl;
  //dom.Particles[3863]->print_history = true;
			
	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		dom.Particles[a]->G		= G;
		dom.Particles[a]->PresEq	= 0;
		dom.Particles[a]->Cs		= Cs;
		dom.Particles[a]->Shepard	= false;
		dom.Particles[a]->Material	= 2;
		//dom.Particles[a]->Et_m = 0.01 * 68.9e9;	//In bilinear this is calculate once, TODO: Change to material definition
		dom.Particles[a]->Et_m = 0.0;	//In bilinear this is calculate once, TODO: Change to material definition
		dom.Particles[a]->Fail		= 1;
		dom.Particles[a]->Sigmay	= Fy;
		dom.Particles[a]->Alpha		= 1.0;
		//dom.Particles[a]->Beta		= 1.0;
		dom.Particles[a]->TI		= 0.3;
		dom.Particles[a]->TIInitDist	= dx;
		double z = dom.Particles[a]->x(2);
    
		// if ( z > L )
			// dom.Particles[a]->ID=3;

    		if ( z > L/2. - dx ){
    			dom.Particles[a]->ID=11;	         
        }
        
	}
  	
	dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
	dom.BC.InOutFlow = 0;
  dom.auto_ts = false;

	dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-5,"test06",1000);
	
	dom.WriteXDMF("ContactTest");
}
MECHSYS_CATCH