#include "Mesh.h"
#include "Domain.h"
#include <iostream>
#include "InteractionAlt.cpp"
#include "SolverKickDrift.cpp"
#include "SolverFraser.cpp"

#define TAU		0.005
#define VMAX	10.0

using namespace SPH;
using namespace std;

std::ofstream of;

double tout;

bool contact = true; //NOT WORKING WITH CONTACT = FALSE

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
  if (contact)
    domi.trimesh[0]->SetVel(Vec3_t(0.0,-vcompress,0.0));


  double dtout = 1.e-4;
  if (domi.getTime()>=tout){
    // cout << "Normal integrated force " <<domi.m_scalar_prop<<endl;
    // cout << "Normal acc sum " << normal_acc_sum<<endl;
    tout += dtout;
    of << domi.getTime()<< ", " << domi.max_disp[2]<<", " << domi.contact_force_sum << ", " << ", " <<domi.ext_forces_work<<", " <<domi.plastic_work << ", " <<domi.accum_cont_heat_cond << ", " << domi.contact_friction_work<<endl;
  }

}


int main(){
	//
	TriMesh mesh;
  mesh.dimension = 2;

	cout << "Creating Mesh" << endl;


	 SPH::Domain	dom;

	dom.Dimension	= 2;
	dom.Nproc	= 12;
	dom.Kernel_Set(Qubic_Spline);
	dom.Scheme	= 1;	//Mod Verlet
	//dom.XSPH	= 0.1; //Very important

		double dx,h,rho,K,G,Cs,Fy;
	double R,L,n;

	R	= 0.15;
	L	= 0.1;
	n	= 30.0;		//in length, radius is same distance

	rho	= 2700.0;
	K	= 6.7549e10;
	G	= 2.5902e10;
	Fy	= 300.e6;
	//dx	= L / (n-1);
	//dx = L/(n-1);
	dx = L/40.;
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

									
	//dom.AddCylinderLength(0, Vec3_t(0.,0.,-L/20.), R, L + 2.*L/20.,  dx/2., rho, h, false); 
  dom.AddBoxLength(0 ,Vec3_t ( 0., 0., 0.0 ), L, L,  0.0 , dx/2.0 ,rho, h, 1 , 0 , false, false );
  
	cout << "Max z plane position: " <<dom.Particles[dom.Particles.Size()-1]->x(1)<<endl;

	double cyl_ymax = dom.Particles[dom.Particles.Size()-1]->x(1) + 1.0800001 * dom.Particles[dom.Particles.Size()-1]->h /*- 1.e-6*/;

	cout << "Plane z " <<cyl_ymax<<endl;
	mesh.AxisPlaneMesh(1,false,Vec3_t(-dx,cyl_ymax, 0.),Vec3_t(0.15,0.0, 0.),50);
	cout << "Plane z" << *mesh.node[0]<<endl;
	
	
	//mesh.AxisPlaneMesh(2,true,Vec3_t(-R-R/10.,-R-R/10.,-L/10.),Vec3_t(R + R/10., R + R/10.,-L/10.),4);
	cout << "Creating Spheres.."<<endl;
	//mesh.v = Vec3_t(0.,0.,);
	mesh.CalcSpheres(); //DONE ONCE

	//////////////////////

	//ALWAYS AFTER SPH PARTICLES
	//TODO: DO THIS INSIDE SOLVER CHECKS
	double hfac = 1.1;	//Used only for Neighbour search radius cutoff
											//Not for any force calc in contact formulation
  cout << "Adding mesh particles ...";
	dom.AddTrimeshParticles(&mesh, hfac, 10); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
  cout << "done."<<endl;
  
	cout << "Done."<<endl;
	dom.ts_nb_inc = 5;
	dom.gradKernelCorr = false; //ATTENTION! USE CFL = 0.7 AND NOT 1.0, IF 1.0 IS USED RESULT DIVERGES
			
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
		dom.Particles[a]->Beta		= 1.0;
		dom.Particles[a]->TI		= 0.3;
		dom.Particles[a]->TIInitDist	= dx;

    double y = dom.Particles[a]->x(1);

		if ( y < 0 ){
			dom.Particles[a]->ID=2;
      dom.Particles[a]->not_write_surf_ID = true;		
		}
	}
	//Contact Penalty and Damping Factors
  if (contact)
    dom.contact = true;
	dom.friction_dyn = 0.0;
	dom.friction_sta = 0.0;
	dom.PFAC = 0.5;
	dom.DFAC = 0.0;
  dom.fric_type = Fr_Dyn;

	
	dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
	dom.BC.InOutFlow = 0;


	of = std::ofstream ("cf.csv", std::ios::out);
    of << "Time, disp, cf, ext_f_wk, plastic_wk, heat_cond, friction_wk"<<endl;
    tout = 0.;  

	//ID 	0 Internal
	//		1	Outer Surface
	//		2,3 //Boundaries
  //dom.auto_ts = false; 
  dom.CFL = 0.7;
  timestep = (0.7*h/(Cs)); //Standard modified Verlet do not accept such step
  //dom.auto_ts=false;

  dom.auto_ts=false;
  //dom.auto_ts_cont = true;
    

  //dom.SolveDiffUpdateFraser(/*tf*/0.0105,/*dt*/timestep,timestep,"test06",1000);
    dom.SolveDiffUpdateFraser(/*tf*/10.0*timestep,/*dt*/timestep,timestep,"test06",1000);
	
	dom.WriteXDMF("ContactTest");
}