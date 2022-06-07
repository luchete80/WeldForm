#include "Mesh.h"
#include "Domain.h"
#include <iostream>
#include "InteractionAlt.cpp"

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
  //cout <<  "Sum "<<domi.m_scalar_prop<<endl;

  
	//TODO: Modify this by relating FEM & AND partciles 
	//domi.trimesh->ApplyConstVel(Vec3_t(0.0,0.0,0.0));
	// domi.trimesh->ApplyConstVel(Vec3_t(0.0,0.0,-vcompress/2.));
  // for (int i = domi.first_fem_particle_idx;i<domi.Particles.Size();i++){
    // domi.Particles[i]->a = Vec3_t(0.0,0.0,0.0);
    // domi.Particles[i]->v = domi.Particles[i]->va = domi.Particles[i]->vb = Vec3_t(0.0,0.0,-vcompress/2.);
  // }
  domi.trimesh[0]->SetVel(Vec3_t(0.0,0.,-vcompress));
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


  double half_plane_length = 0.009;
	int count = 2*0.009/(2.*dx); //Half the density... od original mesh

//double cyl_zmax = L/2. + 4.94e-4; //ORIGINAL
	double cyl_zmax = L/2. + dx*0.545; //If new meshing 
  
  cout << "plane length particle count: "<<count<<endl;
	mesh.AxisPlaneMesh(2,false, Vec3_t(-half_plane_length,-half_plane_length, cyl_zmax),
                              Vec3_t( half_plane_length, half_plane_length, cyl_zmax),count);
	cout << "Plane z" << *mesh.node[0]<<endl;
	
	
	//mesh.AxisPlaneMesh(2,true,Vec3_t(-R-R/10.,-R-R/10.,-L/10.),Vec3_t(R + R/10., R + R/10.,-L/10.),4);
	cout << "Creating Spheres.."<<endl;
	//mesh.v = Vec3_t(0.,0.,);
	mesh.CalcSpheres(); //DONE ONCE
	
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
		dom.Particles[a]->Alpha		= 2.5;
		dom.Particles[a]->Beta		= 2.5;
		dom.Particles[a]->TI		= 0.3;
		dom.Particles[a]->TIInitDist	= dx;
		double z = dom.Particles[a]->x(2);
    
		// if ( z > L )
			// dom.Particles[a]->ID=3;
	}
	//Contact Penalty and Damping Factors
	dom.contact = true;
	dom.friction_dyn = 0.1;
  dom.friction_sta = 0.0;
  dom.fric_type = Fr_Bound;
  
  dom.friction = 0.1;
  
	//dom.friction = 0.0;
	dom.PFAC = 0.3;
	dom.DFAC = 0.2;
  //dom.update_contact_surface = false;
	
	dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
	dom.BC.InOutFlow = 0;

	
	//////////////////////

	//ALWAYS AFTER SPH PARTICLES
	//TODO: DO THIS INSIDE SOLVER CHECKS
	double hfac = 1.1;	//Used only for Neighbour search radius cutoff
											//Not for any force calc in contact formulation
	dom.AddTrimeshParticles(&mesh, hfac, 10); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
	//ID 	0 Internal
	//		1	Outer Surface
	//		2,3 //Boundaries
  dom.auto_ts = false;
	dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-5,"test06",1000);
  
	timestep = (0.4*h/(Cs+VMAX)); //CHANGED WITH VELOCITY
  dom.SolveDiffUpdateKickDrift(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-5,"test06",1000);
	
	dom.WriteXDMF("ContactTest");
}
MECHSYS_CATCH