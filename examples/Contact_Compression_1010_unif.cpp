#include "Mesh.h"
#include "Domain.h"
#include <iostream>
#include "InteractionAlt.cpp"
#include "SolverKickDrift.cpp"
#include "SolverFraser.cpp"

#include "Geometry.cpp"

#define TAU		0.0005
#define VMAX	1.0
#define DX 0.0012

using namespace SPH;
using namespace std;
double tout, dtout;

ofstream ofprop("cont_comp_1010.csv", std::ios::out);
  
int part_per_row;
void UserAcc(SPH::Domain & domi) {
	double vcompress;
		vcompress = VMAX;
	if (domi.getTime() < TAU ) 
		vcompress = VMAX/TAU * domi.getTime();
	else
		vcompress = VMAX;

  double dS = DX * DX;
  domi.m_scalar_prop = 0.;
  double normal_acc_sum=0.;
  for (int i = 0;i<=part_per_row;i++){
    domi.m_scalar_prop += domi.Particles[i]->Sigma (2,2) * dS;
      normal_acc_sum      += domi.Particles[i]->a(2) * domi.Particles[i]->Mass;
  }
  dtout = 1.0e-4;
  if (domi.getTime()>=tout){
    cout << "Normal integrated force " <<domi.m_scalar_prop<<endl;
    cout << "Normal acc sum " << normal_acc_sum<<endl;
    tout += dtout;
    ofprop << domi.max_disp[2]<<", " << domi.contact_force_sum << endl;
  }

  
	// for (int k=0; k<Nproc;k++) 
		// for (size_t a = 0; a < ContPairs[k].Size();a++){
      // if (dot(Particles[a] -> contforce, Particles[a] -> contforce)>0){
        
      // }
    // }
	// #pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	// #ifdef __GNUC__
	// for (size_t i=0; i<domi.Particles.Size(); i++)
	// #else
	// for (int i=0; i<domi.Particles.Size(); i++)
	// #endif
	// {
    // //zcenter, x=R, y=0 particles
		// if (domi.Particles[i]->ID == 2) {
			// domi.Particles[i]->a[1] = domi.Particles[i]->a[2] = 0.;
			// domi.Particles[i]->v[1] = domi.Particles[i]->v[2] = 0.;
			// domi.Particles[i]->va[1] = domi.Particles[i]->va[2] = 0.;
			// domi.Particles[i]->vb[1] = domi.Particles[i]->vb[2] = 0.;

			// //domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		// }    
    // //zcenter, x=0, y=R particles
		// if (domi.Particles[i]->ID == 3) {
			// domi.Particles[i]->a[0] = domi.Particles[i]->a[2] = 0.;
			// domi.Particles[i]->v[0] = domi.Particles[i]->v[2] = 0.;
			// domi.Particles[i]->va[0] = domi.Particles[i]->va[2] = 0.;
			// domi.Particles[i]->vb[0] = domi.Particles[i]->vb[2] = 0.;

			// //domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		// }  
		// if (domi.Particles[i]->ID == 4) {
			// domi.Particles[i]->a = 0.;
			// domi.Particles[i]->v = 0.;
      // domi.Particles[i]->va = 0.;
      // domi.Particles[i]->vb = 0.;

			// //domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		// }        
    // //CENTER TOP AND BOTTOM, IN ORDER TO NOT TO SLIDE IF CONTACT IS UNSTABLE 
		// if (domi.Particles[i]->ID == 5) {
			// domi.Particles[i]->a[0] = domi.Particles[i]->a[1] = 0.;
			// domi.Particles[i]->v[0] = domi.Particles[i]->v[1] = 0.;
			// domi.Particles[i]->va[0] = domi.Particles[i]->va[1] = 0.;
			// domi.Particles[i]->vb[0] = domi.Particles[i]->vb[1] = 0.;

			// //domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		// }     
	// }
  
  domi.trimesh[0]->SetVel(Vec3_t(0.0,0.,-vcompress));
  domi.trimesh[1]->SetVel(Vec3_t(0.0,0., 0.));
}


int main() try{
	//
	TriMesh mesh, mesh2;

	cout << "Creating Mesh" << endl;
  tout = 0.;

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
	//dom.AddCylinderLength(0, Vec3_t(0.,0.,0.), R, L,  dx/2., rho, h, false, ghost); 
  dom.AddCylUniformLength(0, R, L, dx/2,rho, h);
	cout << "Max z plane position: " <<dom.Particles[dom.Particles.Size()-1]->x(2)<<endl;

	double cyl_zmax = dom.Particles[dom.Particles.Size()-1]->x(2) + /*1.005 * */dom.Particles[dom.Particles.Size()-1]->h /*- 1.e-6*/;

	double cyl_zmin = dom.Particles[0]->x(2) - /*1.005 * */dom.Particles[dom.Particles.Size()-1]->h /*- 1.e-6*/;

	mesh.AxisPlaneMesh (2,false,Vec3_t(-1.5*R,-1.5*R, cyl_zmax),Vec3_t(1.5*R,1.5*R, cyl_zmax),20);
  mesh2.AxisPlaneMesh(2,true,Vec3_t(-1.5*R,-1.5*R, cyl_zmin),Vec3_t(1.5*R,1.5*R,cyl_zmin),20);
	
  cout << "Plane z" << *mesh.node[0]<<endl;
	
	
	//mesh.AxisPlaneMesh(2,true,Vec3_t(-R-R/10.,-R-R/10.,-L/10.),Vec3_t(R + R/10., R + R/10.,-L/10.),4);
	cout << "Creating Spheres.."<<endl;
	//mesh.v = Vec3_t(0.,0.,);
	mesh.CalcSpheres(); //DONE ONCE
  mesh2.CalcSpheres();
  
	cout << "Done."<<endl;
	dom.ts_nb_inc = 5;
	dom.gradKernelCorr = false;
	int top, bottom, center;
  top = bottom = center = 0;   
  int center_top = 0;			
    int center_bottom = 0;
  
  double zstart = dom.Particles[0]->x(2);
  bool search = true;
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
		dom.Particles[a]->Alpha		= 1.0;
		dom.Particles[a]->Beta		= 0.0;
		dom.Particles[a]->TI		= 0.3;
		dom.Particles[a]->TIInitDist	= dx;

    double x = dom.Particles[a]->x(0);
    double y = dom.Particles[a]->x(1);
		double z = dom.Particles[a]->x(2);
    
    if (z  > zstart && search ){
      part_per_row = a;
      search = false;
      cout << "Particle per row "<<part_per_row<<endl;
    }
 
    // //If friction is null, the cylinder not slide
    // if ( abs (z - (L/2.-dx)) < dx/2. && (abs(x - R) < 1.5*dx  || abs(x + R) < 1.5*dx ) && abs(y) < 1.1*dx){
      // dom.Particles[a]->ID=2;	  
      // dom.Particles[a]->not_write_surf_ID = true;
      // top++;      
    // } 
    // //x=R, y=0
    // if ( abs (z - (L/2.-dx)) < dx/2. && abs(x) < 1.1*dx && (abs(y-R) < 1.5 *dx || abs(y+R) < 1.5*dx)){
      // dom.Particles[a]->ID=3;	  
      // dom.Particles[a]->not_write_surf_ID = true;
      // bottom++;      
    // } 
    // if ( abs (z - (L/2.-dx)) < dx/2. && abs(x) < dx/2. && abs(y) < dx/2.){
      // dom.Particles[a]->ID=4;	  
      // dom.Particles[a]->not_write_surf_ID = true;
      // center++;      
    // }     

    // if ( z < dx/2. && abs(x) < dx/2. && abs(y) < dx/2.){
      // dom.Particles[a]->ID=5;	  
      // dom.Particles[a]->not_write_surf_ID = true;
      // center_bottom++;      
    // } 

    // if ( z > (L - 1.5*dx) && abs(x) < dx/2. && abs(y) < dx/2.){
      // dom.Particles[a]->ID=5;	  
      // dom.Particles[a]->not_write_surf_ID = true;
      // center_top++;      
    // }   
    
    // if ( z > (L - 1.5*dx) && abs(x*x+y*y) > R-2.*dx){
      // //dom.Particles[a]->ID=5;	  
      // //dom.Particles[a]->not_write_surf_ID = true;
      // center_top++;      
    // }   
	}
  cout << top<< " Side 1 particles, "<<bottom << " side 2 particles, "<<center << " center particles" <<endl; 
  cout << "Center PERIMETER: " <<center_top <<endl;
  cout << "Center Bottom: " <<center_bottom <<endl;  
	//Contact Penalty and Damping Factors
	dom.contact = true;
	dom.friction_dyn = 0.2;
  dom.friction_sta = 0.2;
  dom.fric_type = Fr_Dyn;
 
	dom.PFAC = 0.6;
	dom.DFAC = 0.0;

	//ALWAYS AFTER SPH PARTICLES
	//TODO: DO THIS INSIDE SOLVER CHECKS
	double hfac = 1.1;	//Used only for Neighbour search radius cutoff
											//Not for any force calc in contact formulation
	dom.AddTrimeshParticles(&mesh, hfac, 10); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
  dom.AddTrimeshParticles(&mesh2, hfac, 11); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
	//ID 	0 Internal
	//		1	Outer Surface
	//		2,3 //Boundaries
  dom.auto_ts     = false;
  //dom.auto_ts_acc = true;
//  	dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-5,"test06",1000);
  
	// timestep = (0.4*h/(Cs+VMAX)); //CHANGED WITH VELOCITY
  // dom.SolveDiffUpdateKickDrift(/*tf*/0.105,/*dt*/timestep,/*dtOut*/1.e-4,"test06",1000);
  dom.auto_ts=true; 
  dom.CFL = 0.7;
  timestep = (0.7*h/(Cs+VMAX)); //CHANGED WITH VELOCITY
  
  //dom.SolveDiffUpdateLeapfrog(/*tf*/0.105,/*dt*/timestep,/*dtOut*/1.e-4,"test06",1000);  
	dom.SolveDiffUpdateFraser(/*tf*/0.10105,/*dt*/timestep,/*dtOut*/1.e-4,"test06",1000);
  //dom.SolveDiffUpdateFraser(/*tf*/0.01205,/*dt*/timestep,/*dtOut*/timestep,"test06",1000);  
	
  dom.WriteXDMF("ContactTest");
}
MECHSYS_CATCH