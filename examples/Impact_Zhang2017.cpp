
#include <iostream>
#include "Domain.h"
#include "Input.h"

#include "Input.h"
#include "InteractionAlt.cpp"
#include "Mesh.h"

#include "SolverFraser.cpp"
#include "Geometry.cpp"
#include "SolverKickDrift.cpp"
#include "SolverLeapfrog.cpp"
#define TAU		0.005
#define VMAX	190.0

using namespace SPH;
using namespace std;


//////////////////////////////////
///// EXAMPLE FROM 
 
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
		if (domi.Particles[i]->ID == 10) // "FEM", fictitious SPH PARTICLES FROM TRIMESH
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,-vcompress);
			domi.Particles[i]->va		= Vec3_t(0.0,0.0,-vcompress);
//			domi.Particles[i]->vb		= Vec3_t(0.0,0.0,-vcompress);
//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
		if (domi.Particles[i]->ID == 2)
		{
			// domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			// domi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);
			// domi.Particles[i]->vb		= Vec3_t(0.0,0.0,0.0);
			//domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
	}
	
	//TODO: Modify this by relating FEM & AND partciles 
	//domi.trimesh->ApplyConstVel(Vec3_t(0.0,0.0,0.0));
	domi.trimesh[0]->ApplyConstVel(Vec3_t(0.0,0.0,-vcompress));
}


int main(){
	//
	TriMesh mesh;

	cout << "Creating Mesh" << endl;


	 SPH::Domain	dom;

	dom.Dimension	= 3;
	dom.Nproc	= 4;
	dom.Kernel_Set(Qubic_Spline);
	dom.Scheme	= 1;	//Mod Verlet
	//dom.XSPH	= 0.1; //Very important

	double dx,h,rho,E,nu,K,G,Cs,Fy;
	double R,L,n;
	//J-Cook material
	//A 175
	//B 380
	//C 0.0015
	//n 0.34
	//m 1
	// Material A (MPa) B (MPa) C n m Troom (K) Tmelt (K
// OFHC copper 90 292 0.0025 0. 31 1.09 273 1356
// Armco iron 175 380 0.06 0.32 0.55 273 1811

	R	= 0.0032;
	L	= 0.0324;
	n	= 60.0;		//in length, radius is same distance

	rho	= 8930.0;
	
	E = 117.0e9;
	nu = 0.35;
	K = E / ( 3.*(1.-2*nu) );
	G = E / (2.* (1.+nu));

	Elastic_ el(E,nu);
	//Hollomon(const double eps0_, const double &k_, const double &m_):
	//Hollomon mat(el,Fy/E,1220.e6,0.195);
	double eps_0 = 1.0;
	//Zhang 2017
	//Smoothed particle hydrodynamics with kernel gradient correction for
	//modeling high velocity impact in two- and three-dimensional spaces
	//Material A (MPa) B (MPa) C n m Troom (K) Tmelt (K)
	//OFHC copper 90 292 0.0025 0.31 1.09 273 1356
	//Armco iron 175 380 0.06 0.32 0.55 273 1811

	double A,B,C,n_,m,T_m,T_t;
	A = 90.e6; B = 292.0e6; C = 0.0025;
	m = 1.09;  n_ = 0.31; eps_0 = 1.0;
	T_m = 1356.; T_t = 273.;

	JohnsonCook mat(el, A,B,C,
                      m,n_,eps_0,
                      T_m, T_t);	
											
	
	Fy	= 400.e6;
	dx	= L / (n-1);
	//dx = 0.001;
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

									
	dom.AddCylinderLength(0, Vec3_t(0.,0.,-L/10.), R, L + 2.*L/10.,  dx/2., rho, h, false); 
	cout << "Max z plane position: " <<dom.Particles[dom.Particles.Size()-1]->x(2)<<endl;

	double cyl_zmax = dom.Particles[dom.Particles.Size()-1]->x(2) + 1.000001 * dom.Particles[dom.Particles.Size()-1]->h /*- 1.e-6*/;

	
	mesh.AxisPlaneMesh(2,false,Vec3_t(-0.01,-0.01, cyl_zmax),Vec3_t(0.01,0.01, cyl_zmax),40);
	cout << "Plane z" << *mesh.node[0]<<endl;
	
	
	//mesh.AxisPlaneMesh(2,true,Vec3_t(-R-R/10.,-R-R/10.,-L/10.),Vec3_t(R + R/10., R + R/10.,-L/10.),4);
	cout << "Creating Spheres.."<<endl;
	//mesh.v = Vec3_t(0.,0.,);
	mesh.CalcSpheres(); //DONE ONCE
	
	double T_h = 0.;	//Homologous or room temp
	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		dom.Particles[a]->mat = &mat;
		dom.Particles[a]->T = 273.;
		dom.Particles[a]-> Material_model = JOHNSON_COOK;
		//Need to calculate Initial Yield Stress
		//CalcYieldStress(const double &strain, const double &strain_rate, const double &temp)
		dom.Particles[a]->Sigmay = dom.Particles[a]->mat->CalcYieldStress(0.,0., 273.0);	//For  Hollomon

		dom.Particles[a]->G		= G;
		dom.Particles[a]->PresEq	= 0;
		dom.Particles[a]->Cs		= Cs;
		dom.Particles[a]->Shepard	= false;
		dom.Particles[a]->Material	= 2;
		//dom.Particles[a]->Et_m = 0.01 * 68.9e9;	//In bilinear this is calculate once, TODO: Change to material definition
		dom.Particles[a]->Et_m = 0.0;	//In bilinear this is calculate once, TODO: Change to material definition
		dom.Particles[a]->Fail		= 1;
		//dom.Particles[a]->Sigmay	= Fy;
		dom.Particles[a]->Alpha		= 1.5;
		dom.Particles[a]->Beta		= 1.5;
		dom.Particles[a]->TI		= 0.3;
		dom.Particles[a]->TIInitDist	= dx;
		double z = dom.Particles[a]->x(2);
		if ( z < 0 ){
			dom.Particles[a]->ID=2;
			dom.Particles[a]->IsFree=false;
			dom.Particles[a]->NoSlip=true;			
		
		}
		// if ( z > L )
			// dom.Particles[a]->ID=3;
	}
	//Contact Penalty and Damping Factors
	dom.contact = true;
	dom.friction = 0.0;
	dom.PFAC = 1.0;
	dom.DFAC = 0.2;
	dom.update_contact_surface = false;
	
	dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
	dom.BC.InOutFlow = 0;

	
	//////////////////////

	//ALWAYS AFTER SPH PARTICLES
	//TODO: DO THIS INSIDE SOLVER CHECKS
	double hfac = 1.1;	//Used only for Neighbour search radius cutoff
											//Not for any force calc in contact formulation
	dom.AddTrimeshParticles(&mesh, hfac, 10); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
  
  timestep = (0.6*h/(Cs+VMAX)); 
  dom.CFL = 0.6;
  //timestep = 2.5e-6;
  dom.auto_ts = true;
	//ID 	0 Internal
	//		1	Outer Surface
	//		2,3 //Boundaries
	//dom.Solve(/*tf*/40.e-6,/*dt*/timestep,/*dtOut*/1.e-6,"test06",1000);
	//dom.Solve(/*tf*/60.01e-6,/*dt*/timestep,/*dtOut*/1.0e-6,"test06",999);
  //dom.SolveDiffUpdateFraser(/*tf*/60.0e-6,/*dt*/timestep,/*dtOut*/1.e-6,"test06",1000);
  dom.SolveDiffUpdateLeapFrog(/*tf*/60.0e-6,/*dt*/timestep,/*dtOut*/1.e-6,"test06",1000);
	//dom.ThermalStructSolve(/*tf*/60.01e-6,/*dt*/timestep,/*dtOut*/1.0e-6,"test06",999);
	
	dom.WriteXDMF("ContactTest");
}