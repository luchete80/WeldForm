#include "Mesh.h"
#include "Domain.h"
#include <iostream>
#include "InteractionAlt.cpp"
#include "SolverFraser.cpp"
#include "SolverKickDrift.cpp"

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
    //In fraser algorithm acceleration should be fixed
		if (domi.Particles[i]->ID == 2) {
			domi.Particles[i]->v			= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->a			= Vec3_t(0.0,0.0,0.0);
    }
    else if (domi.Particles[i]->ID == 3) {
			domi.Particles[i]->v			= Vec3_t(0.0,0.0,-vcompress);
			domi.Particles[i]->a			= Vec3_t(0.0,0.0,0.0);
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

	bool ghost = false;								
	dom.AddCylinderLength(0, Vec3_t(0.,0.,0.), R, L,  dx/2., rho, h, false, ghost); 
	cout << "Max z plane position: " <<dom.Particles[dom.Particles.Size()-1]->x(2)<<endl;

	dom.gradKernelCorr = false;

	int top_part = 0;
  int bottom_part = 0;
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
    
		if ( z > L-dx){
      top_part++;
			dom.Particles[a]->ID=3;
    }
    else if ( z < dx){
      bottom_part++;
			dom.Particles[a]->ID=2;
    }
	}
  cout << "top_part: "<<top_part<<endl;
  cout << "bottom_part: "<<bottom_part<<endl;
  
//  	dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-5,"test06",1000);
  dom.auto_ts=false; 
  dom.CFL = 0.6;
  timestep = (0.6*h/(Cs+VMAX)); //CHANGED WITH VELOCITY
  
  //dom.SolveDiffUpdateKickDrift(/*tf*/0.105,/*dt*/timestep,/*dtOut*/1.e-4,"test06",1000);
	//dom.SolveDiffUpdateLeapfrog(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-4 ,"test06",1000);
  dom.SolveDiffUpdateFraser(/*tf*/0.01205,/*dt*/timestep,/*dtOut*/1.e-4,"test06",1000);
	
  dom.WriteXDMF("ContactTest");
}
=======
#include "Mesh.h"
#include "Domain.h"
#include <iostream>
#include "InteractionAlt.cpp"
#include "SolverFraser.cpp"
#include "SolverKickDrift.cpp"

#define VMAX	1.0

using namespace SPH;
using namespace std;

#define DX 0.0012
double tout, dtout;

ofstream ofprop("cont_comp_1010.csv", std::ios::out);

void UserAcc(SPH::Domain & domi) {
	double vcompress;
  double dS = DX * DX;
  double sigma_sum = 0.;
		vcompress = VMAX;
    double normal_acc_sum=0.;
	//#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	//#ifdef __GNUC__
	//for (size_t i=0; i<domi.Particles.Size(); i++)
	//#else
	for (int i=0; i<domi.Particles.Size(); i++)
	//#endif
	
	{
    //In fraser algorithm acceleration should be fixed
		if (domi.Particles[i]->ID == 2) {
			domi.Particles[i]->v			= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->a			= Vec3_t(0.0,0.0,0.0);
      sigma_sum += domi.Particles[i]->Sigma (2,2) * dS;
    }
    else if (domi.Particles[i]->ID == 3) {
			domi.Particles[i]->v			= Vec3_t(0.0,0.0,-vcompress);
			domi.Particles[i]->a			= Vec3_t(0.0,0.0,0.0);
    }
  }
  
  dtout = 1.0e-4;
  if (domi.getTime()>=tout){
    tout += dtout;
    cout << "Time: "<<  domi.getTime()<<"Sigma sum "<<sigma_sum<<endl;
    ofprop << domi.max_disp[2]<<", " <<  domi.m_scalar_prop << endl;
  }
}


int main() try{
	//
	TriMesh mesh;
  tout = 0.;
  
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

	dx = DX;  //Tenth of radius
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

	bool ghost = false;								
	dom.AddCylinderLength(0, Vec3_t(0.,0.,0.), R, L,  dx/2., rho, h, false, ghost); 
	cout << "Max z plane position: " <<dom.Particles[dom.Particles.Size()-1]->x(2)<<endl;

	dom.gradKernelCorr = true;

	int top_part = 0;
  int bottom_part = 0;
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
    
		if ( z > L-dx){
      top_part++;
			dom.Particles[a]->ID=3;
    }
    else if ( z < dx){
      bottom_part++;
			dom.Particles[a]->ID=2;
    }
	}
  cout << "top_part: "<<top_part<<endl;
  cout << "bottom_part: "<<bottom_part<<endl;
  
//  	dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-5,"test06",1000);
  dom.auto_ts=false; 
  dom.CFL = 0.6;
  timestep = (0.6*h/(Cs+VMAX)); //CHANGED WITH VELOCITY
  
  //dom.SolveDiffUpdateKickDrift(/*tf*/0.105,/*dt*/timestep,/*dtOut*/1.e-4,"test06",1000);
	//dom.SolveDiffUpdateLeapfrog(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-4 ,"test06",1000);
  dom.SolveDiffUpdateFraser(/*tf*/0.01205,/*dt*/timestep,/*dtOut*/1.e-4,"test06",1000);
	
  dom.WriteXDMF("ContactTest");
}
MECHSYS_CATCH