#include "Mesh.h"
#include "Domain.h"
#include "InteractionAlt.cpp"

#include <iostream>

#define TAU		0.005
#define VMAX	0.1 

using namespace SPH;
using namespace std;

std::ofstream of;

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
  int center = 0;
  int side = 0;
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	#ifdef __GNUC__
	for (size_t i=0; i<domi.Particles.Size(); i++)
	#else
	for (int i=0; i<domi.Particles.Size(); i++)
	#endif	
	{
    //Center particles
		if (domi.Particles[i]->ID == 2) {
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->vb		= Vec3_t(0.0,0.0,0.0);
      center++;
			//domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
    //zcenter, x=R, y=0 particles
		if (domi.Particles[i]->ID == 3) {
			domi.Particles[i]->a[1] = domi.Particles[i]->a[2] = 0.;
			domi.Particles[i]->v[1] = domi.Particles[i]->v[2] = 0.;
			domi.Particles[i]->va[1] = domi.Particles[i]->va[2] = 0.;
			domi.Particles[i]->vb[1] = domi.Particles[i]->vb[2] = 0.;
      side++;
			//domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}    
  }
  //cout << "center particles, "<<center<< ", side particles"<<side<<endl;
  //To avoid spurious slide
  int id = domi.Particles.Size()-1;
  domi.Particles[id]->v [0]  =   domi.Particles[id]->v[1]  = 0.;
  domi.Particles[id]->va[0] =   domi.Particles[id]->va[1] = 0.;
  domi.Particles[id]->vb[0] =   domi.Particles[id]->vb[1] = 0.;

  domi.Particles[0]->v[0]  =   domi.Particles[0]->v[1]  = 0.;
  domi.Particles[0]->va[0] =   domi.Particles[0]->va[1] = 0.;
  domi.Particles[0]->vb[0] =   domi.Particles[0]->vb[1] = 0.;
  
  domi.trimesh[0]->SetVel(Vec3_t(0.0,0.,-vcompress));
  domi.trimesh[1]->SetVel(Vec3_t(0.0,0., 0.));
  of << domi.getTime() << ", "<<domi.Particles[7451]->contforce(2)<<endl;
  //cout << "Position "<<domi.Particles[7451]->x<<endl;
}


int main() try{
	//
	TriMesh mesh,mesh2;

	cout << "Creating Mesh" << endl;


	 SPH::Domain	dom;

	dom.Dimension	= 3;
	dom.Nproc	= 4;
	dom.Kernel_Set(Qubic_Spline);
	dom.Scheme	= 2;	//Mod Verlet
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

  double Et = (1400.e6-Fy)/0.5;
  double Ep =  E*Et/(E-Et);
  //Ep = 0.;
	//dx	= L / (n-1);
	//dx = L/(n-1);
	dx = 0.0008;  //Tenth of radius
	h	= dx*1.2; //Very important
	Cs	= sqrt(K/rho);

	double timestep;


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
  //dom.AddFixedMassScaling(1000.);
	cout << "Max z plane position: " <<dom.Particles[dom.Particles.Size()-1]->x(2)<<endl;


  double half_plane_length = 0.012;
	int count = 2*0.009/(2.*dx); //Half the density... od original mesh

//double cyl_zmax = L/2. + 4.94e-4; //ORIGINAL
	double cyl_zmax = L - 0.1*h; //If new meshing 
  
  cout << "plane length particle count: "<<count<<endl;
	mesh.AxisPlaneMesh(2,false, Vec3_t(-half_plane_length,-half_plane_length, cyl_zmax),
                              Vec3_t( half_plane_length, half_plane_length, cyl_zmax),count);
	cout << "Plane z" << *mesh.node[0]<<endl;

	mesh2.AxisPlaneMesh(2,true, Vec3_t(-half_plane_length,-half_plane_length, -h),
                              Vec3_t( half_plane_length, half_plane_length, -h),count);  

	//mesh.AxisPlaneMesh(2,true,Vec3_t(-R-R/10.,-R-R/10.,-L/10.),Vec3_t(R + R/10., R + R/10.,-L/10.),4);
	cout << "Creating Spheres.."<<endl;
	//mesh.v = Vec3_t(0.,0.,);
	mesh.CalcSpheres(); //DONE ONCE
  mesh2.CalcSpheres();
	
	dom.ts_nb_inc = 5;
	dom.gradKernelCorr = true;
  
  dom.Particles[3863]->print_history = true;
  
  cout << "Particle 3589 coords xyz "<<dom.Particles[3863]->x<<endl;
  //dom.Particles[3863]->print_history = true;
	int top, bottom;
  top = bottom =0;  
	for (size_t a=0; a<dom.Particles.Size(); a++)
	{
		dom.Particles[a]->G		= G;
		dom.Particles[a]->PresEq	= 0;
		dom.Particles[a]->Cs		= Cs;
		dom.Particles[a]->Shepard	= false;
		dom.Particles[a]->Material	= 2;
    dom.Particles[a]->Ep		= Ep;
		//dom.Particles[a]->Et_m = 0.01 * 68.9e9;	//In bilinear this is calculate once, TODO: Change to material definition
		dom.Particles[a]->Et_m = 0.0;	//In bilinear this is calculate once, TODO: Change to material definition
		dom.Particles[a]->Fail		= 1;
		dom.Particles[a]->Sigmay	= Fy;
		dom.Particles[a]->Alpha		= 2.5;
		dom.Particles[a]->Beta		= 2.5;
		dom.Particles[a]->TI		= 0.3;
		dom.Particles[a]->TIInitDist	= dx;
    double x = dom.Particles[a]->x(0);
    double y = dom.Particles[a]->x(1);
		double z = dom.Particles[a]->x(2);
    
		// if ( z > L )
			// dom.Particles[a]->ID=3;
    
    //If friction is null, the cylinder not slide
    if ( abs (z - (L/2.)) < dx/2. && abs(x) < dx/2. && abs(y) < dx/2.){
      dom.Particles[a]->ID=2;	  
      dom.Particles[a]->not_write_surf_ID = true;
      top++;      
    } 
    //x=R, y=0
    if ( abs (z - (L/2.)) < dx/2. && abs(x - R) < 1.5*dx && abs(y) < 1.1*dx){
      dom.Particles[a]->ID=3;	  
      dom.Particles[a]->not_write_surf_ID = true;
      bottom++;      
    } 
    
	}

	//type definition to shorten coding
	std::ostringstream oss;
	of = std::ofstream ("cf.csv", std::ios::out);
  of << "Time, cfpart7451"<<endl;
	//of << oss.str();

  
	dom.contact = true;
	//dom.friction = 0.15;
	dom.friction_dyn = 0.15;
	dom.friction_sta = 0.15;
	dom.PFAC = 0.8;
	dom.DFAC = 0.2;
  dom.fric_type = Fr_Bound;

  cout << top<< " Top particles, "<<bottom << " bottom particles"<<endl; 
  	
	dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
	dom.BC.InOutFlow = 0;
  dom.auto_ts = false;
  //timestep=1.e-8;

	double hfac = 1.1;	//Used only for Neighbour search radius cutoff
  
	dom.AddTrimeshParticles(&mesh, hfac, 10); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
  dom.AddTrimeshParticles(&mesh2, hfac, 11); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
    
  timestep = (0.4*h/(Cs+VMAX)); //Standard modified Verlet do not accept such step
	
  for (size_t a=0; a<dom.Particles.Size(); a++) {
    if (dom.Particles[a]->ID == 2)
      cout << endl<<"ID 2"<<endl;
  }
  //dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-5,"test06",1000);
  dom.SolveDiffUpdateKickDrift(/*tf*/0.105,/*dt*/timestep,/*dtOut*/1.e-5,"test06",10000);
	//dom.SolveDiffUpdateModEuler(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-5,"test06",1000);
  //dom.SolveDiffUpdateModVerlet(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-5,"test06",1000);
	dom.WriteXDMF("ContactTest");
}
MECHSYS_CATCH