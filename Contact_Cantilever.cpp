#include "Domain.h"

#define TAU		0.005
#define UMAX	0.0024

using namespace SPH;
using namespace std;
//Generated curved surface
void GenerateMesh(TriMesh *mesh, const double &r, const double &xcenter, const double dens){
	double ang=-Util::PI/6;
	double anginc = (Util::PI/3.)/dens;
	int n = 0;
	for (int j=0;j<dens+1;j++){
		ang=-Util::PI/6;
		for (int i=0;i<dens+1;i++){
			(*mesh->node[n])[0] = xcenter + r*sin(ang);
			//cout << "orig y"<< (*mesh->node[n])[1] <<endl;
			(*mesh->node[n])[1] -= r*cos(ang);
			//cout << "curr y"<< (*mesh->node[n])[1] <<endl;
			//dom.Particles[i]->x(1) = ycenter - r*cos();
			ang += anginc; //xhange angle in x axis
			n++;
		}	
		//cout <<"z change"<<endl;
	}
	mesh->CalcCentroids();
	mesh->CalcNormals();
}

void UserAcc(SPH::Domain & domi)
{
	double vtraction;

	if (domi.getTime() < TAU ) 
		vtraction = UMAX/TAU ;
	else
		vtraction = 0.0;
	//cout << "vtraction"<<vtraction<<endl;
	
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	#ifdef __GNUC__
	for (size_t i=0; i<domi.Particles.Size(); i++)
	#else
	for (int i=0; i<domi.Particles.Size(); i++)
	#endif
	
	{
		if (domi.Particles[i]->ID == 10)
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,-vtraction,0.0);
			domi.Particles[i]->va		= Vec3_t(0.0,-vtraction,0.0);
//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
		// if (domi.Particles[i]->ID == 2)
		// {
			// domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			// domi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);
			// domi.Particles[i]->vb		= Vec3_t(0.0,0.0,0.0);
// //			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		// }
	}
	domi.trimesh->ApplyConstVel(Vec3_t(0.0,-vtraction,0.0));
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
	//     	dom.XSPH	= 0.5; //Very important

	double dx,h,rho,K,G,Cs,Fy;
	double H,Lx,Ly,Lz;

	Lx = 0.1;
	Ly = 0.024;
	Lz = 0.012;	

	double E  = 70.e9;
	double nu = 0.33;

	rho	= 2700.0;
	K= E / ( 3.*(1.-2*nu) );
	G= E / (2.* (1.+nu));
	Fy	= 1000.e10;

	dx = 0.002;
	h	= dx*1.2; //Very important
	Cs	= sqrt(K/rho);

	double timestep;
	timestep = (0.2*h/(Cs));

	// timestep = 1.0e-9;
	// dom.auto_ts = false;
	//timestep = 5.e-7;

	cout<<"t  = "<<timestep<<endl;
	cout<<"Cs = "<<Cs<<endl;
	cout<<"K  = "<<K<<endl;
	cout<<"G  = "<<G<<endl;
	cout<<"Fy = "<<Fy<<endl;
	dom.GeneralAfter = & UserAcc;
	dom.DomMax(0) = Lx;
	dom.DomMin(0) = -Lx;
									
	dom.AddBoxLength(0 ,Vec3_t ( -Lx/2.0 -Lx/20.0, -Ly/2.0 , -Lz/2.0 ), 
					Lx , Ly /*+dx*/,  Lz /*+dx*/, 
					dx/2.0 ,rho, h, 1 , 0 , false, false );

	double ymax = dom.Particles[dom.Particles.Size()-1]->x(1) + dom.Particles[dom.Particles.Size()-1]->h - 1.2e-4;
	cout << "y max "<< ymax << endl;
	double xmax = dom.Particles[dom.Particles.Size()-1]->x(0);
	TriMesh mesh;
	int mechdens = 8;
	mesh.AxisPlaneMesh(1,true,Vec3_t(xmax-0.01,ymax + Ly,-0.01),Vec3_t(xmax+0.01,ymax + Ly,0.01),mechdens);
	//Change plane coordinates to curved mesh
	GenerateMesh(&mesh, Ly, xmax,mechdens); //This should be done before calc spheres

	mesh.CalcSpheres(); //DONE ONCE
	cout << "Normal, pplane, radius"<<endl;
	for (int e=0;e<mesh.element.Size();e++){
		//cout << mesh.element[e]->normal<<endl;
		//cout << mesh.element[e]->pplane<<endl;
		//cout << mesh.element[e]->radius<<endl;
	}
	
		//ALWAYS AFTER SPH PARTICLES
	//TODO: DO THIS INSIDE SOLVER CHECKS
	double hfac = 1.1;	//Used only for Neighbour search radius cutoff
											//Not for any force calc in contact formulation
	dom.AddTrimeshParticles(mesh, hfac, 10); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
	
	dom.contact = true;	
	dom.friction = 0.0;
	dom.PFAC = 1.0;
	dom.DFAC = 0.0;
	
	//cout << "Plane z" << *mesh.node[0]<<endl;
	
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
    		dom.Particles[a]->TI		= 0.3;
    		dom.Particles[a]->TIInitDist	= dx;
    		double x = dom.Particles[a]->x(0);
			double y = dom.Particles[a]->x(1);
			double z = dom.Particles[a]->x(2);
			
    		if ( x <= -Lx/2. ){
    			dom.Particles[a]->ID=2;
    			dom.Particles[a]->IsFree=false;
    			dom.Particles[a]->NoSlip=true;
			}
    		// if ( y >= (Ly/2. - dx ) && x >= (Lx/2. -Lx/40.) )
    			// dom.Particles[a]->ID=3;
    	}
		dom.WriteXDMF("maz");
//		dom.m_kernel = SPH::iKernel(dom.Dimension,h);	


    	dom.Solve(/*tf*/0.0101,/*dt*/timestep,/*dtOut*/0.00001,"test06",999);
        return 0;
}
MECHSYS_CATCH
