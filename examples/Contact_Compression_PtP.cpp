#include "Mesh.h"
#include "Domain.h"
#include <iostream>

#define TAU		0.005
#define VMAX	10.0

using namespace SPH;
using namespace std;

void UserAcc(SPH::Domain & domi) {
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
		//TODO: Modify this by relating FEM & AND partciles 
		if (domi.Particles[i]->ID == 11) // "FEM", fictitious SPH PARTICLES FROM TRIMESH
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,-vcompress);
			domi.Particles[i]->va		= Vec3_t(0.0,0.0,-vcompress);
			domi.Particles[i]->vb		= Vec3_t(0.0,0.0,-vcompress);
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
	domi.trimesh->ApplyConstVel(Vec3_t(0.0,0.0,-vcompress));
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

									
	dom.AddCylinderLength(0, Vec3_t(0.,0.,-L/10.), R, L + 2.*L/10. ,  dx/2., rho, h, false); 
	
	dom.first_fem_particle_idx = dom.Particles.Size(); //For Surface Calculation
	
	
	double rho2 = 3.*rho;
	// inline void Domain::AddBoxLength(int tag, Vec3_t const & V, double Lx, double Ly, double Lz, 
									// double r, double Density, double h, int type, int rotation, bool random, bool Fixed)
									
	dom.AddBoxLength(10, Vec3_t (-1.2*R,-1.2*R,L + L/10. - dx/2. /*-dx/4.*/), 2.4*R, 2.4*R, L/6., 
									dx/2., rho2, h, 1 , 0 , false, false );
									
	cout << "Max z plane position: " <<dom.Particles[dom.Particles.Size()-1]->x(2)<<endl;

	double cyl_zmax = dom.Particles[dom.Particles.Size()-1]->x(2) + 1.000001 * dom.Particles[dom.Particles.Size()-1]->h /*- 1.e-6*/;

	
	//mesh.AxisPlaneMesh(2,false,Vec3_t(-0.5,-0.5, cyl_zmax),Vec3_t(0.5,0.5, cyl_zmax),40);
	//cout << "Plane z" << *mesh.node[0]<<endl;
	
	
	//mesh.AxisPlaneMesh(2,true,Vec3_t(-R-R/10.,-R-R/10.,-L/10.),Vec3_t(R + R/10., R + R/10.,-L/10.),4);
	//cout << "Creating Spheres.."<<endl;
	//mesh.v = Vec3_t(0.,0.,);
	//mesh.CalcSpheres(); //DONE ONCE
	
	dom.ts_nb_inc = 1;

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
		if ( z < 0 ){
			dom.Particles[a]->ID=2;
			dom.Particles[a]->IsFree=false;
			dom.Particles[a]->NoSlip=true;			
		
		}
		
		if (z > (L + L/10.+L/10.)){
			dom.Particles[a]->ID=11;			
		}
		// if ( z > L )
			// dom.Particles[a]->ID=3;
	}

	//Change Second Body
	double E = 210.0e10;	//
	double nu = 0.3;
	K= E / ( 3.*(1.-2*nu) );
	G= E / (2.* (1.+nu));
	Cs	= sqrt(K/rho);
	
	for (size_t a=dom.first_fem_particle_idx; a<dom.Particles.Size(); a++) {
		dom.Particles[a]->Cs		= Cs;
		dom.Particles[a]->G		= 10.*G;
		dom.Particles[a]->Sigmay	= 1000.*Fy;
	}
		
		
	//Contact Penalty and Damping Factors
	dom.contact = true;
	dom.friction = 0.15;
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
	//dom.AddTrimeshParticles(mesh, hfac, 10); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
	//ID 	0 Internal
	//		1	Outer Surface
	//		2,3 //Boundaries
	dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-5,"test06",1000);
	
	dom.WriteXDMF("ContactTest");
}