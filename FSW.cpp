#include "Mesh.h"
#include "Domain.h"
#include "NastranReader.h"


#define VFAC			10.0
#define VAVA			35.			//mm/min
#define WROT 			1200.0 	//rpm

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


	}
	Vec3_t omega(0.,WROT*M_PI/30.*VFAC);
  domi.trimesh->ApplyConstVel(Vec3_t(0.0,0.0,-VAVA * VFAC));
	domi.trimesh->RotateAxisVel(omega, domi.getStepSize());
	
  for (int i = domi.first_fem_particle_idx;i<domi.Particles.Size();i++){
    domi.Particles[i]->a = Vec3_t(0.0,0.0,0.0);
    domi.Particles[i]->v = domi.Particles[i]->va = domi.Particles[i]->vb = Vec3_t(0.0,0.0, - VAVA * VFAC);
  }

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

		H	= 0.005;
		L	= 0.05;
    
		n	= 30.0;		//in length, radius is same distance

	rho	= 2700.0;
	K	= 6.7549e10;
	G	= 2.5902e10;
	Fy	= 300.e6;
	//dx	= L / (n-1);
	//dx = L/(n-1);
	dx = 0.00085;
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


		// inline void Domain::AddCylinderLength(int tag, Vec3_t const & V, double Rxy, double Lz, 
									// double r, double Density, double h, bool Fixed) {
		bool symlength = true; //Also Symmetric on z axis
		bool Fixed = false;
	//Cylinder Slice of 90 degree and half length
	// THIS DOES NOT HAVE Z INITIAL POSITION SINCE IT IS ZERO OR -LZ ACCORGIND TO LAST ARGUMENT
	// void AddDoubleSymCylinderLength(int tag, double Rxy, double Lz, 
									// double r, double Density, double h, bool Fixed, bool symlength = false);
  
  dom.AddBoxLength(0 ,Vec3_t ( -L/2.0-L/20.0 , -H -H/6 , -L/2.0-L/20.0 ), L + L/10.0 + dx/10.0 , H ,  L + L/10. , dx/2.0 ,rho, h, 1 , 0 , false, false );

  SPH::NastranReader reader("Tool.nas");
  
  SPH::TriMesh mesh(reader);
  
  //double cyl_zmax = L/2. + 4.94e-4; //ORIGINAL
  double cyl_zmax = L/2. + dx*0.6 -1.e-3; //If new meshing  

	cout << "Creating contact mesh.."<<endl;
	
	cout << "Plane z" << *mesh.node[0]<<endl;
  
  cout << "Mesh node size "<<mesh.node.Size()<<endl;
	
	
	//mesh.AxisPlaneMesh(2,true,Vec3_t(-R-R/10.,-R-R/10.,-L/10.),Vec3_t(R + R/10., R + R/10.,-L/10.),4);
	cout << "Creating Spheres.."<<endl;
	//mesh.v = Vec3_t(0.,0.,);
	mesh.CalcSpheres(); //DONE ONCE
	double hfac = 1.1;	//Used only for Neighbour search radius cutoff
											//Not for any force calc in contact formulation
	dom.AddTrimeshParticles(mesh, hfac, 10); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
    
  
	dom.ts_nb_inc = 5;
	dom.gradKernelCorr = false;
			
	cout << "Particle count: "<<dom.Particles.Size()<<endl;

		for (size_t a=0; a<dom.Particles.Size(); a++)
		{
			dom.Particles[a]->G		= G;
			dom.Particles[a]->PresEq	= 0;
			dom.Particles[a]->Cs		= Cs;
			dom.Particles[a]->Shepard	= false;
			dom.Particles[a]->Material	= 2;
			dom.Particles[a]->Fail		= 1;
			dom.Particles[a]->Sigmay	= Fy;
			dom.Particles[a]->Alpha		= 1.0;
			//dom.Particles[a]->Beta		= 1.0;
			dom.Particles[a]->TI		= 0.3;
			dom.Particles[a]->TIInitDist	= dx;
			
			double x = dom.Particles[a]->x(0);
			double y = dom.Particles[a]->x(1);
			double z = dom.Particles[a]->x(2);
			
			
			//BOTTOM PLANE
			// if ( z < dx  && z > -dx/2. ){
				// dom.Particles[a]->ID=3;
				// dom.Particles[a]->not_write_surf_ID = true;
			// }
			
			
			// if ( x < dx  && x > -dx/2. && z < L/2. - dx)
				// dom.Particles[a]->ID=1;
			// if ( y < dx  && y > -dx/2. && z < L/2. - dx)
				// dom.Particles[a]->ID=2; 
				
			//x,y, central symmetry
			// if ( y < dx  && y > -dx/2. && x < dx  && x > -dx/2. && z < L/2. - dx)
				// dom.Particles[a]->ID=4;  
			
			// if ( y < dx  && y > -dx/2. && z < dx  && z > -dx/2. ) //yz -5
				// dom.Particles[a]->ID=5;           
			// if (  x < dx  && x > -dx/2. && z < dx  && z > -dx/2. ) //xz - 6
				// dom.Particles[a]->ID=6;         
			
			//First one captures 4 particle, second one
			//if ( y < dx  && y > -dx && x < dx  && x > -dx && z < dx  && z > -dx/2. ) //xz - 7
			// if ( y < dx  && y > -dx && x < dx  && x > -dx/2. && z < dx  && z > -dx/2. ) //xz - 7
				// dom.Particles[a]->ID=7;   

			//TOP
			// if ( y < dx  && y > -dx/2. && z > L/2. - dx ) //yz -5
				// dom.Particles[a]->ID=8;           
			// if (  x < dx  && x > -dx/2. && z > L/2. - dx ) //xz - 6
				// dom.Particles[a]->ID=9;         
			// if ( y < dx  && y > -dx/2. && x < dx  && x > -dx/2. && z > L/2. - dx ) //xyz - 7
				// dom.Particles[a]->ID=10;         
		}
		
    // dom.Particles[0]->IsFree=false;
    // dom.Particles[0]->NoSlip=true;			
	//Contact Penalty and Damping Factors
  dom.fric_type = Fr_Dyn;
	dom.contact = true;
	//dom.friction = 0.15;
	dom.friction_dyn = 0.1;
  dom.friction_sta = 0.0;
	dom.PFAC = 0.8;
	dom.DFAC = 0.2;
	dom.update_contact_surface = false;

  dom.WriteXDMF("maz");
  dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
  dom.BC.InOutFlow = 0;

  dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/0.0001,"test06",999);
  
  return 0;
}
MECHSYS_CATCH
