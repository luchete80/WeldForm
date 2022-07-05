#include "Mesh.h"
#include "Domain.h"
#include <iostream>
#include "InteractionAlt.cpp"

#define TAU		0.005
#define VMAX	20.0

using namespace SPH;
using namespace std;

std::ofstream of;
#define TWO_SURF 1

void UserAcc(SPH::Domain & domi) {
	double vcompress = VMAX;

	// if (domi.getTime() < TAU ) 
		// vcompress = VMAX/TAU * domi.getTime();
	// else
		// vcompress = VMAX;
	//cout << "time: "<< domi.getTime() << "V compress "<< vcompress <<endl;
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)

	#ifdef __GNUC__
	for (size_t i=0; i<domi.Particles.Size(); i++)
	#else
	for (int i=0; i<domi.Particles.Size(); i++)
	#endif
	
	{
		//TODO: Modify this by relating FEM & AND partciles 
		// if (domi.Particles[i]->ID == 10) // "FEM", fictitious SPH PARTICLES FROM TRIMESH
		// {
			// domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			// domi.Particles[i]->v		= Vec3_t(0.0,0.0,-vcompress);
			// domi.Particles[i]->va		= Vec3_t(0.0,0.0,-vcompress);
// //			domi.Particles[i]->vb		= Vec3_t(0.0,0.0,-vcompress);
// //			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		// }
		if (domi.Particles[i]->ID == 2)
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->vb		= Vec3_t(0.0,0.0,0.0);
			//domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
	}
	
	//TODO: Modify this by relating FEM & AND partciles 
	//domi.trimesh->ApplyConstVel(Vec3_t(0.0,0.0,0.0));
	//domi.trimesh->ApplyConstVel(Vec3_t(0.0,0.0,-vcompress));
  domi.trimesh[0]->SetVel(Vec3_t(0.0,0.,-vcompress));
  //#ifdef TWO_SURF
  domi.trimesh[1]->SetVel(Vec3_t(0.0,0.,0));
  //#endif
  of << domi.getTime() << ", "<<domi.Particles[12419]->contforce(2)<< ", " << domi.Particles[12419]->v(2)<<endl;
}


int main() try{
	//
	TriMesh mesh,mesh2;

	cout << "Creating Mesh" << endl;


	 SPH::Domain	dom;

	dom.Dimension	= 3;
	dom.Nproc	= 4;
	dom.Kernel_Set(Qubic_Spline);
	dom.Scheme	= 1;	//Mod Verlet
	//dom.XSPH	= 0.1; //Very important

		double dx,h,rho,K,G,Cs,Fy;
	double R,L,n;

	R	= 0.127;
	L	= 0.3;
	n	= 30.0;		//in length, radius is same distance

	rho	= 7850.0;
  double E = 200.e9;
  double nu = 0.3;
  K= E / ( 3.*(1.-2*nu) );
  G= E / (2.* (1.+nu));

	Elastic_ el(E,nu);
  Fy	= 260.e6;
  Hollomon mat(el,Fy,7.1568e8,0.22);
  
	//dx	= L / (n-1);
	//dx = L/(n-1);
	dx = 0.012;
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

									
	dom.AddCylinderLength(0, Vec3_t(0.,0.,0.), R, L,  dx/2., rho, h, false); 
	cout << "Max z plane position: " <<dom.Particles[dom.Particles.Size()-1]->x(2)<<endl;

	double cyl_zmax = dom.Particles[dom.Particles.Size()-1]->x(2) + 1.000001 * dom.Particles[dom.Particles.Size()-1]->h /*- 1.e-6*/;

	cout << "Plane z " <<cyl_zmax<<endl;
	mesh.AxisPlaneMesh(2,false,Vec3_t(-0.15,-0.15, cyl_zmax),Vec3_t(0.15,0.15, cyl_zmax),20);
	cout << "Plane z" << *mesh.node[0]<<endl;
	
  //#ifdef TWO_SURF
  mesh2.AxisPlaneMesh(2,false,Vec3_t(-0.15,-0.15, -h),Vec3_t(0.15,0.15, -h),20);
	//#endif
	//mesh.AxisPlaneMesh(2,true,Vec3_t(-R-R/10.,-R-R/10.,-L/10.),Vec3_t(R + R/10., R + R/10.,-L/10.),4);
	cout << "Creating Spheres.."<<endl;
	//mesh.v = Vec3_t(0.,0.,);
	mesh.CalcSpheres(); //DONE ONCE
  
  //#ifdef TWO_SURF
  mesh2.CalcSpheres();
  //#endif
	
  cout << "Done."<<endl;
	dom.ts_nb_inc = 5;
	dom.gradKernelCorr = true;
			
	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		dom.Particles[a]->G		= G;
		dom.Particles[a]->PresEq	= 0;
		dom.Particles[a]->Cs		= Cs;
		dom.Particles[a]->Shepard	= false;
		dom.Particles[a]->Material	= 2;
    dom.Particles[a]-> Material_model = HOLLOMON;
    dom.Particles[a]->mat = &mat;
    
		//dom.Particles[a]->Et_m = 0.01 * 68.9e9;	//In bilinear this is calculate once, TODO: Change to material definition
		dom.Particles[a]->Et_m = 0.0;	//In bilinear this is calculate once, TODO: Change to material definition
		dom.Particles[a]->Fail		= 1;
		dom.Particles[a]->Sigmay	= Fy;
		dom.Particles[a]->Alpha		= 1.0;
		//dom.Particles[a]->Beta		= 1.0;
		dom.Particles[a]->TI		= 0.3;
		dom.Particles[a]->TIInitDist	= dx;
		double z = dom.Particles[a]->x(2);
    #ifndef TWO_SURF
		if ( z < dx ){
			dom.Particles[a]->ID=2;	
      dom.Particles[a]->not_write_surf_ID = true;		
		}
    #endif
		// if ( z > L )
			// dom.Particles[a]->ID=3;
	}
	//Contact Penalty and Damping Factors
	dom.contact = true;
	dom.friction_dyn = 0.15;
	dom.friction_sta = 0.2;
	dom.PFAC = 0.8;
	dom.DFAC = 0.0;
  dom.fric_type = Fr_Bound;

	
	dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
	dom.BC.InOutFlow = 0;


	of = std::ofstream ("cf.csv", std::ios::out);
  of << "Time, cf, vypart"<<endl;

	//////////////////////

	//ALWAYS AFTER SPH PARTICLES
	//TODO: DO THIS INSIDE SOLVER CHECKS
	double hfac = 1.1;	//Used only for Neighbour search radius cutoff
											//Not for any force calc in contact formulation
  cout << "Adding mesh particles ...";
	dom.AddTrimeshParticles(&mesh, hfac, 10); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
  //#ifdef TWO_SURF
  //dom.AddTrimeshParticles(&mesh2, hfac, 11); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
  //#endif
  cout << "done."<<endl;
    

	//ID 	0 Internal
	//		1	Outer Surface
	//		2,3 //Boundaries
  //dom.auto_ts = false;
  timestep = (0.4*h/(Cs)); //Standard modified Verlet do not accept such step
  //dom.auto_ts=false;

  dom.auto_ts=true;
    
	//dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-5,"test06",1000);
  //THIS DOES NOT WORK WITH FIXED PARTICLES
  dom.SolveDiffUpdateKickDrift(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-4 ,"test06",1000);
	
	dom.WriteXDMF("ContactTest");
}

MECHSYS_CATCH