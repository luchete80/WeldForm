#include "Mesh.h"
#include "Domain.h"
#include <iostream>
#include "InteractionAlt.cpp"
#include "SolverKickDrift.cpp"
#include "SolverFraser.cpp"

#define VMAX	0.1

using namespace SPH;
using namespace std;

void UserAcc(SPH::Domain & domi) {
	double vcompress;
		vcompress = VMAX;
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	#ifdef __GNUC__
	for (size_t i=0; i<domi.Particles.Size(); i++)
	#else
	for (int i=0; i<domi.Particles.Size(); i++)
	#endif
	
	{
		if (domi.Particles[i]->ID == 2) 
			domi.Particles[i]->v			= Vec3_t(0.0,0.0,-vcompress);
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
	dom.Scheme	= 0;	//Mod Verlet
	//dom.XSPH	= 0.1; //Very important

		double dx,h,rho,K,G,Cs,Fy;
	double R,L;

	R	= 0.0127;
	L	= 0.030;

	rho	= 7850.0;
  double E = 200.e9;
  double nu = 0.3;
  K= E / ( 3.*(1.-2*nu) );
  G= E / (2.* (1.+nu));

	dx = 0.0012;  //Tenth of radius
	h	= dx*1.2; //Very important
	Cs	= sqrt(K/rho);

	Elastic_ el(E,nu);
  Fy	= 260.e6;
  Hollomon mat(el,Fy,7.1568e8,0.22);
  
	double timestep;
	timestep = (0.1*h/(Cs+VMAX)); //CHANGED WITH VELOCITY

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

	dom.gradKernelCorr = false;

	int top_part = 0;
	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		dom.Particles[a]->G		= G;
		dom.Particles[a]->PresEq	= 0;
		dom.Particles[a]->Cs		= Cs;
		dom.Particles[a]->Shepard	= false;
    dom.Particles[a]-> Material_model = HOLLOMON;
    dom.Particles[a]->mat = &mat;
        
		dom.Particles[a]->Fail		= 1;
		dom.Particles[a]->Sigmay	= Fy;
		dom.Particles[a]->Alpha		= 2.5;
		dom.Particles[a]->Beta		= 2.5;
		dom.Particles[a]->TI		= 0.3;
		dom.Particles[a]->TIInitDist	= dx;
		double z = dom.Particles[a]->x(2);
    
		if ( z > L/2. - 2.0 * dx){
      top_part++;
			dom.Particles[a]->ID=2;
    }
	}
  cout << "top_part: "<<top_part<<endl;

  dom.auto_ts = false;
//  	dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-5,"test06",1000);
  dom.CFL = 0.7; //For auto ts
	timestep = (0.7*h/(Cs+VMAX)); //CHANGED WITH VELOCITY
  //dom.SolveDiffUpdateKickDrift(/*tf*/0.105,/*dt*/timestep,/*dtOut*/1.e-4,"test06",1000);
  //dom.SolveDiffUpdateLeapfrog(/*tf*/0.105,/*dt*/timestep,/*dtOut*/1.e-4,"test06",1000);  
  dom.SolveDiffUpdateFraser(/*tf*/0.01205,/*dt*/timestep,/*dtOut*/1.e-4,"test06",1000);
    
	dom.WriteXDMF("ContactTest");
}
MECHSYS_CATCH