#include "Mesh.h"
#include "Domain.h"
#include "NastranReader.h"
#include "InteractionAlt.cpp"
#include "SolverKickDrift.cpp"
#include "SolverFraser.cpp"

#define VFAC			15.0
#define VAVA			5.833e-4		//35 mm/min
#define WROT 			1200.0 	    //rpm
#define TOOLRAD   0.0062
#define SUPPRAD   0.01

void UserAcc(SPH::Domain & domi) {
	double vcompress;

	//cout << "time: "<< domi.getTime() << "V compress "<< vcompress <<endl;
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)

	#ifdef __GNUC__
	for (size_t i=0; i<domi.Particles.Size(); i++)
	#else
	for (int i=0; i<domi.Particles.Size(); i++)
	#endif
	
	{
    //Vertical Constraint
		if (domi.Particles[i]->ID == 2) {
			//domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
      domi.Particles[i]->a(1)		  = 0.;
			domi.Particles[i]->va(1)	  = 0.;
			domi.Particles[i]->v(1)		  = 0.;
			domi.Particles[i]->vb(1)	  = 0.;
      domi.Particles[i]->VXSPH(1) = 0.;
			//domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
    
    if (domi.Particles[i]->ID == 3) {
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);      
    }

	}
	//Mesh is updated automatically
  // for (int i = domi.first_fem_particle_idx;i<domi.Particles.Size();i++){
    // domi.Particles[i]->a = Vec3_t(0.0,0.0,0.0);
    // //THIS SHOULD BE FROM BARICENTER AND AUTOMATIC!
    // //domi.Particles[i]->v = domi.Particles[i]->va = domi.Particles[i]->vb = Vec3_t(0.0,- VAVA * VFAC,0.);
  // }
  //domi.trimesh[0]->SetVel(Vec3_t(0.0,0.,-vcompress));
}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
		 SPH::Domain	dom;

		dom.Dimension	= 3;
		dom.Nproc	= 4;
		dom.Kernel_Set(Qubic_Spline);
		dom.Scheme	= 1;	//Mod Verlet
		//dom.XSPH	= 0.1; //Very important

    double dx,h,rho,K,G,Cs,Fy;
		double H,L,n;

		H	= 0.001;
		L	= 0.05;
    
		n	= 30.0;		//in length, radius is same distance

	rho	= 2700.0;
	// K	= 6.7549e10;
	// G	= 2.5902e10;
	double E  = 70.e9;
	double nu = 0.3;

    
	K= E / ( 3.*(1.-2*nu) );
	G= E / (2.* (1.+nu));

	Fy	= 350.e6;  
	//Fy	= 300.e6;
	//dx	= L / (n-1);
	//dx = L/(n-1);
	dx = 0.00025;
	h	= dx*1.2; //Very important
	Cs	= sqrt(K/rho);

	double timestep;
	timestep = (0.2*h/(Cs));

	Elastic_ el(E,nu);
///// MATERIAL CONSTANTS EXAMPLE FROM
///// Zhang_2017 (Aluminium) ?
  double A,B,C,n_,m,T_m,T_t,eps_0;
  A = 276.e6; B = 255.0e6; C = 0.0015;
  m = 1.0;  n_ = 0.3; eps_0 = 1.0;
  T_m = 775.; T_t = 273.;
			
  // 𝐴
  // 𝐽𝐶 276.0 MPa
  // 𝐵
  // 𝐽𝐶 255.0 MPa
  // 𝑛𝐽𝐶 0.3 -
  // 𝑚𝐽𝐶 1.0

	JohnsonCook mat(el, A,B,C,
                      m,n_,eps_0,
                      T_m, T_t);	
                      
		//timestep = 2.5e-6;

	cout<<"t  = "<<timestep<<endl;
	cout<<"Cs = "<<Cs<<endl;
	cout<<"K  = "<<K<<endl;
	cout<<"G  = "<<G<<endl;
	cout<<"Fy = "<<Fy<<endl;
	dom.GeneralAfter = & UserAcc;
	dom.DomMax(0) = L;
	dom.DomMin(0) = -L;


		// inline void Domain::AddCylinderLength(int tag, Vec3_t const & V, double Rxy, double Lz, 
									// double r, double Density, double h, bool Fixed) {
		bool symlength = true; //Also Symmetric on z axis
		bool Fixed = false;
	//Cylinder Slice of 90 degree and half length
	// THIS DOES NOT HAVE Z INITIAL POSITION SINCE IT IS ZERO OR -LZ ACCORGIND TO LAST ARGUMENT
	// void AddDoubleSymCylinderLength(int tag, double Rxy, double Lz, 
									// double r, double Density, double h, bool Fixed, bool symlength = false);
  
  double ybottom = -H - dx/1.45;  
  
  double ytop = ybottom + H ; 
  
  dom.AddBoxLength(0 ,Vec3_t ( -L/2.0-L/20.0 , ybottom, -L/2.0-L/20.0 ), L + L/10.0 + dx/10.0 , H ,  L + L/10., dx/2.0 ,rho, h, 1 , 0 , false, false );

  SPH::NastranReader reader("Tool.nas");
  
  SPH::TriMesh mesh(reader);
  
  //double cyl_zmax = L/2. + 4.94e-4; //ORIGINAL
  double cyl_zmax = L/2. + 5.0 * dx/*-1.e-3*/; //If new meshing  

	cout << "Creating contact mesh.."<<endl;
	
	cout << "Plane z" << *mesh.node[0]<<endl;
  
  cout << "Mesh node size "<<mesh.node.Size()<<endl;
	
	
	//mesh.AxisPlaneMesh(2,true,Vec3_t(-R-R/10.,-R-R/10.,-L/10.),Vec3_t(R + R/10., R + R/10.,-L/10.),4);
	cout << "Creating Spheres.."<<endl;
	//mesh.v = Vec3_t(0.,0.,);
	mesh.CalcSpheres(); //DONE ONCE
	double hfac = 1.1;	//Used only for Neighbour search radius cutoff
											//Not for any force calc in contact formulation
	dom.AddTrimeshParticles(&mesh, hfac, 10); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
    
  
	dom.ts_nb_inc = 5;
	dom.gradKernelCorr = false;
			
	cout << "Particle count: "<<dom.Particles.Size()<<endl;
  int bottom_particles = 0;
  int top_particles = 0;
  int side_particles  =0;
		for (size_t a=0; a<dom.Particles.Size(); a++)
		{
      dom.Particles[a]-> Material_model = JOHNSON_COOK/*HOLLOMON*/;
      dom.Particles[a]->mat = &mat;
        
			dom.Particles[a]->G		= G;
			dom.Particles[a]->PresEq	= 0;
			dom.Particles[a]->Cs		= Cs;
			dom.Particles[a]->Shepard	= false;
			dom.Particles[a]->Material	= 2;
			dom.Particles[a]->Fail		= 1;
			dom.Particles[a]->Sigmay	= Fy;
			dom.Particles[a]->Alpha		= 1.0;
			//dom.Particles[a]->Beta		= 2.0;
			dom.Particles[a]->TI		= 0.3;
			dom.Particles[a]->TIInitDist	= dx;
      
      dom.Particles[a]->k_T			  =	130.*VFAC;  //[W/(m.K)]
			dom.Particles[a]->cp_T			=	960.;       //[J/(kg.K)]
			// dom.Particles[a]->h_conv		= 100.0; //W/m2-K
			// dom.Particles[a]->T_inf 		= 500.;
			dom.Particles[a]->T				  = 20.0;			
			
			double x = dom.Particles[a]->x(0);
			double y = dom.Particles[a]->x(1);
			double z = dom.Particles[a]->x(2);
			
      double r = sqrt (x*x+z*z);      
      if (/*r < TOOLRAD && */y < (ybottom +dx ) ){
        dom.Particles[a]->ID=2; //ID 1 is free surface  
        dom.Particles[a]->not_write_surf_ID = true;
        bottom_particles++;
        
        if (r < TOOLRAD){
          dom.Particles[a]->ID = 4; //ID 1 is free surface  
          dom.Particles[a]->not_write_surf_ID = true;
          dom.Particles[a]->Thermal_BC 	= TH_BC_CONVECTION;
          dom.Particles[a]->h_conv		= 200.0 * VFAC; //W/m2-K
          dom.Particles[a]->T_inf 		= 20.;
        }
      }

      if (r > SUPPRAD && y > ( ytop - dx ) ){
        dom.Particles[a]->ID=2; //ID 1 is free surface  
        dom.Particles[a]->not_write_surf_ID = true;
        top_particles++;
      }
			
			
			//SIDES
			if ( z < -L/2. + dx || z > L/2. - dx){
				dom.Particles[a]->ID=3;
				dom.Particles[a]->not_write_surf_ID = true;
   			dom.Particles[a]->IsFree=false;
        side_particles++;
			}
			else if ( x < -L/2. + 2.*dx || x > L/2. - 2.*dx){
				dom.Particles[a]->ID=3;
				dom.Particles[a]->not_write_surf_ID = true;
        dom.Particles[a]->IsFree=false;
        side_particles++;
			}			
     
		}
    
  cout << "Bottom particles: " << bottom_particles << endl;
  cout << "Top particles: " << top_particles << endl;
  cout << "Side particles: " << side_particles<<endl;
		
    // dom.Particles[0]->IsFree=false;
    // dom.Particles[0]->NoSlip=true;			
	//Contact Penalty and Damping Factors
  dom.fric_type = Fr_Dyn;
	dom.contact = true;
	//dom.friction = 0.15;
	dom.friction_dyn = 0.2;
  dom.friction_sta = 0.2;
	dom.PFAC = 0.6;
	dom.DFAC = 0.0;
	dom.update_contact_surface = false;

  dom.WriteXDMF("maz");
  dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
  dom.BC.InOutFlow = 0;
  
  // SET TOOL BOUNDARY CONDITIONS
  //dom.trimesh[0]->SetRotAxisVel(Vec3_t(0.,WROT*M_PI/30.*VFAC,0.));  //axis rotation m_w
  dom.trimesh[0]->SetVel(Vec3_t(0.0,-VAVA * VFAC,0.));              //translation, m_v


  dom.auto_ts = false;        //AUTO TS FAILS IN THIS PROBLEM (ISSUE)
  dom.thermal_solver = true;
  dom.cont_heat_gen = true;
   
  timestep = (0.7*h/(Cs)); //Standard modified Verlet do not accept such step
  //dom.auto_ts=false;

  dom.auto_ts=true;
  
  //dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/400* timestep,"test06",999);
  dom.SolveDiffUpdateFraser(/*tf*/0.1,/*dt*/timestep,/*dtOut*/timestep  ,"test06",1000);
    
  return 0;
}
MECHSYS_CATCH
