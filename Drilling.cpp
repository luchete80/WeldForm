#include "Mesh.h"
#include "Domain.h"
#include "NastranReader.h"
#include "InteractionAlt.cpp"
#include "SolverKickDrift.cpp"
#include "SolverFraser.cpp"
#include "Geometry.cpp" //CylinderUniformLength

bool bottom_contact = false;

#define VFAC			10.0 
//#define VAVA			5.833e-4		//35 mm/min
#define VAVA			14.0e-3		//35 mm/min
#define WROT 			3600.0 	    //rpm

#define ANG_CORTE   10.0  //degrees
#define Z_TIP       0.001 //
#define TAU         0.005 //

#define   TTOT       2.0e-2
#define  DTOUT       1.0e-4 //WITHOUT VFAC
//drill bit points to z positive 

double tout, dtout;

ofstream ofprop("drill_force.csv", std::ios::out);

bool cf_cte = true;

void UserAcc(SPH::Domain & domi) {
	double vava, wrot;

 
	//cout << "time: "<< domi.getTime() << "V compress "<< vcompress <<endl;
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
  #ifdef __GNUC__
  for (size_t i=0; i<domi.Particles.Size(); i++)
  #else
  for (int i=0; i<domi.Particles.Size(); i++)
  #endif	
  {
    
    //BOTTOM
    if (domi.Particles[i]->ID == 2) {
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

  
	if (domi.getTime() <(TAU /VFAC)){ 
		vava  = VAVA/TAU * domi.getTime()         * (VFAC*VFAC);
    wrot  = WROT/TAU *domi.getTime() *M_PI/30.* (VFAC*VFAC);
	} else {
    vava  = VAVA * VFAC;
    wrot  = WROT*M_PI/30.*VFAC;
  }

  if (domi.getTime()>=tout){
    tout += DTOUT/VFAC;
    ofprop << domi.max_disp[2]<<", " << domi.contact_force_sum << endl;
    cout << "V ava: "<<vava<<endl;
  }
  
  /////// SET TOOL BOUNDARY CONDITIONS
  domi.trimesh[0]->SetRotAxisVel(Vec3_t(0.,0.0,wrot));  //axis rotation m_w
  domi.trimesh[0]->SetVel(Vec3_t(0.0,0.0,vava));              //translation, m_v
  
}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
  

    
    
  SPH::Domain	dom;

  dom.Dimension	= 3;
  dom.Nproc	= 32;
  dom.Kernel_Set(Qubic_Spline);
  dom.Scheme	= 1;	//Mod Verlet
  //dom.XSPH	= 0.1; //Very important
  if (argc > 1) {
    string input=argv[1];	
    int proc = atoi(input.c_str());
    dom.Nproc	= proc;
  }
  
  cout << "Threads "<<dom.Nproc	<<endl;
  double dx,h,rho,K,G,Cs,Fy;
  double L,R, n;

  L	= 0.002;
  R = 0.004;

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
	dx = 0.00008;  
	h	= dx*1.2; //Very important
	Cs	= sqrt(K/rho);

	double timestep;
	timestep = (0.2*h/(Cs));

	Elastic_ el(E,nu);
///// MATERIAL CONSTANTS EXAMPLE FROM
///// Zhang_2017 (Aluminium) ?
  double A,B,C,n_,m,T_m,T_t,eps_0;
  A = 352.0e6; B = 440.0e6; C = 0.0083;
  m = 1.0;  n_ = 0.42; eps_0 = 1.0;
  T_m = 502.; T_t = 0.;
			
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
  
  //double ybottom = -H - 1.2 * dx;  /////LARGE PIN, Original Tool
  double ybottom = -L - 0.7 * dx;  /////SMALL PIN, New Tool, IF dx = 0.0005; //If dx = 0.0006 then 
  
  //double ybottom = -H +0.3*dx; //if dx = 0.0002
  
  //TODO: make gap adjustment automatic
  
  double ytop = ybottom + L ; 
  
  // void Domain::AddCylUniformLength(int tag, double Rxy, double Lz, 
																				// double r, double Density, double h) {
  //Rxy LZ
  dom.AddCylUniformLength(0., R, L, dx/2.0, rho, h);
  cout << "Created. "<<endl;
  ///////////////////// DELETE MATERIAL 
  double angle = ANG_CORTE * M_PI / 180.0;
  
  for (size_t i=0; i<dom.Particles.Size(); i++) {
    double r = sqrt(pow(dom.Particles[i]->x(0), 2.0) + pow(dom.Particles[i]->x(1),2.0));
    if (r < 0.00375){

      if (dom.Particles[i]->x(2) < Z_TIP - r*tan(angle)){
        dom.Particles[i]->ID = 2;
      }
    } else {//radius
      if (dom.Particles[i]->x(2) < Z_TIP - 0.00375 * tan(angle)){
        dom.Particles[i]->ID = 2;
      }
    }
  }
  // // // MASS RECALC??
  // // //PRINT TOTMASS
  dom.DelParticles (2); //
  
  cout << "New mesh size: "<<dom.Particles.Size()<<endl;
  dom.solid_part_count = dom.Particles.Size();
  
  
  SPH::NastranReader reader("mecha.nas");
  // //SPH::NastranReader reader("Tool.nas");
  
  SPH::TriMesh mesh(reader/*,true*/); //flip normals
  // SPH::TriMesh mesh2;//Only if bottom contact
  
  // //double cyl_zmax = L/2. + 4.94e-4; //ORIGINAL
  // double cyl_zmax = L/2. + 5.0 * dx/*-1.e-3*/; //If new meshing  
  cout <<"Moving mesh"<<endl;
	// cout << "Creating Spheres.."<<endl;
	// //mesh.v = Vec3_t(0.,0.,);
	mesh.CalcSpheres(); //DONE ONCE, BEFORE ANY MOVE!
  mesh.Move(Vec3_t(0.0,0.0,Z_TIP - 1.2 *h)); //Original Z is zero : CRASH
  	// for (int n=0;n<mesh.node.Size();n++){
		// *mesh.node[n] += Vec3_t(0.0,0.0,2.0*Z_TIP - h);
	// } 
  
	cout << "Creating contact mesh.."<<endl;
	
	cout << "Plane z" << *mesh.node[0]<<endl;
  
  cout << "Mesh node size "<<mesh.node.Size()<<endl;
 
	

	double hfac = 1.1;	//Used only for Neighbour search radius cutoff
											//Not for any force calc in contact formulation
	dom.AddTrimeshParticles(&mesh, hfac, 10); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){

 
	dom.ts_nb_inc = 5;
	dom.gradKernelCorr = false;
  dom.nonlock_sum = true;
			
	cout << "Particle count: "<<dom.Particles.Size()<<endl;
  int bottom_particles = 0;
  int top_particles = 0;
  int side_particles  =0;
		for (size_t a=0; a<dom.solid_part_count; a++)
		{
      dom.Particles[a]-> Material_model = JOHNSON_COOK/*HOLLOMON*/;
      dom.Particles[a]->mat = &mat;
        
			dom.Particles[a]->G		= G;
			dom.Particles[a]->PresEq	= 0;
			dom.Particles[a]->Cs		= Cs;
			dom.Particles[a]->Shepard	= false;
			dom.Particles[a]->Material	= 2;
			dom.Particles[a]->Fail		= 1;

      dom.Particles[a]->Sigmay	= mat.CalcYieldStress(0.0,0.0,0.);    
			
      dom.Particles[a]->Alpha		= 1.5;
			dom.Particles[a]->Beta		= 0.6;
			dom.Particles[a]->TI		= 0.4;
			dom.Particles[a]->TIInitDist	= dx;
      
      dom.Particles[a]->k_T			  =	121.*VFAC;  //[W/(m.K)]
			dom.Particles[a]->cp_T			=	875.;       //[J/(kg.K)]
			// dom.Particles[a]->h_conv		= 100.0; //W/m2-K
			// dom.Particles[a]->T_inf 		= 500.;
			dom.Particles[a]->T				  = 20.0;			
			
			double x = dom.Particles[a]->x(0);
			double y = dom.Particles[a]->x(1);
			double z = dom.Particles[a]->x(2);
			
      double r = sqrt (x*x+y*y);      
      if (r > R - dx ){
        dom.Particles[a]->ID=2; //ID 1 is free surface  
        dom.Particles[a]->not_write_surf_ID = true;
        bottom_particles++;
        

        // dom.Particles[a]->Thermal_BC 	= TH_BC_CONVECTION;
        // dom.Particles[a]->h_conv		= 1000.0 * VFAC; //W/m2-K
        // dom.Particles[a]->T_inf 		= 20.;
      }
      	
			         
		}
    
  cout << "Bottom particles: " << bottom_particles << endl;

    // dom.Particles[0]->IsFree=false;
    // dom.Particles[0]->NoSlip=true;			
	//Contact Penalty and Damping Factors
  dom.fric_type = Fr_Dyn;
	dom.contact = true;
	//dom.friction = 0.15;
  
  if (cf_cte) {
    dom.friction_dyn = 0.2;
    dom.friction_sta = 0.2;
  } else {
    dom.friction_function = Linear;
    dom.friction_b = 0.3;
    dom.friction_m = (0.1-0.3)/500.;
  }
	dom.PFAC = 0.6;
	dom.DFAC = 0.0;
	dom.update_contact_surface = false;
  
  cout << "Writing output "<<endl;
  dom.WriteXDMF("maz");
  cout << "Done. "<<endl;

  dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
  dom.BC.InOutFlow = 0;
  
  // // // SET TOOL BOUNDARY CONDITIONS
  // // dom.trimesh[0]->SetRotAxisVel(Vec3_t(0.,0.0,WROT*M_PI/30.*VFAC));  //axis rotation m_w
  // // dom.trimesh[0]->SetVel(Vec3_t(0.0,0.0,VAVA * VFAC));              //translation, m_v
  
  dom.auto_ts = false;        //AUTO TS FAILS IN THIS PROBLEM (ISSUE)
  dom.thermal_solver = true;
  dom.cont_heat_gen = true;
  
  //dom.cont_heat_cond  = true;
  //dom.contact_hc      = 1000.*VFAC; 
   
  timestep = (0.7*h/(Cs)); //Standard modified Verlet do not accept such step

  //dom.SolveDiffUpdateFraser(/*tf*/TTOT/VFAC,/*dt*/timestep,/*dtOut*/DTOUT/VFAC  ,"test06",1000);
  dom.SolveDiffUpdateFraser(/*tf*/0.01,/*dt*/timestep,/*dtOut*/timestep  ,"test06",1000);
    
  return 0;
}
MECHSYS_CATCH
