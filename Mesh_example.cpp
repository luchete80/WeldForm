#include "Mesh.h"
#include "Domain.h"
#include <iostream>

using namespace SPH;
using namespace std;

void UserAcc(SPH::Domain & domi) {
	// double vcompress;

	// if (domi.getTime() < TAU ) 
		// vcompress = VMAX/TAU * domi.getTime();
	// else
		// vcompress = VMAX;
	// //cout << "time: "<< domi.getTime() << "V compress "<< vcompress <<endl;
	// #pragma omp parallel for schedule (static) num_threads(domi.Nproc)

	// #ifdef __GNUC__
	// for (size_t i=0; i<domi.Particles.Size(); i++)
	// #else
	// for (int i=0; i<domi.Particles.Size(); i++)
	// #endif
	
	// {
		// if (domi.Particles[i]->ID == 3)
		// {
			// domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			// domi.Particles[i]->v		= Vec3_t(0.0,0.0,-vcompress);
			// domi.Particles[i]->va		= Vec3_t(0.0,0.0,-vcompress);
			// //domi.Particles[i]->vb		= Vec3_t(0.0,0.0,-vcompress);
// //			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		// }
		// if (domi.Particles[i]->ID == 2)
		// {
			// // domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			// // domi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);
			// // domi.Particles[i]->vb		= Vec3_t(0.0,0.0,0.0);
			// //domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		// }
	// }
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
	dx = 0.018;
	h	= dx*1.1; //Very important
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

	mesh.AxisPlaneMesh(2,true,Vec3_t(-0.5,-0.5,L + L/10.-dx/2.),Vec3_t(0.5,0.5, L + L/10.-dx/2.),30);
	
	
	for (int v=0;v<mesh.node.Size();v++){
			
	}
	
	//mesh.AxisPlaneMesh(2,true,Vec3_t(-R-R/10.,-R-R/10.,-L/10.),Vec3_t(R + R/10., R + R/10.,-L/10.),4);
	cout << "Creating Spheres.."<<endl;
	mesh.CalcSpheres();
	
									
	dom.AddCylinderLength(1, Vec3_t(0.,0.,-L/10.), R, L + 2.*L/10.,  dx/2., rho, h, false); 

	//ALWAYS AFTER SPH PARTICLES
	dom.AddTrimeshParticles(mesh, 1.1, 10); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
		
		
	cout << "Generating Nb List" <<endl; 
	//dom.InitialChecks();
	dom.CellInitiate();
	dom.ListGenerate();
	// dom.PrintInput(TheFileKey);
	// dom.TimestepCheck();
	// dom.WholeVelocity();
	cout << "Main Nb Search"<<endl;
	dom.MainNeighbourSearch_Ext();	//MUST DO CELL INITIATE FIRST
	dom.SaveNeighbourData();				//Necesary to calulate surface! Using Particle->Nb (count), could be included in search
	dom.CalculateSurface(2);				//After Nb search	
	cout << "Contact Nb Search" <<endl;
	
	cout << "FSM Pairs" <<endl;
	for (int k = 0 ; k< dom.Nproc; k++)
		cout << dom.FSMPairs[k].Size() <<endl;
	
	dom.ContactNbSearch();
	dom.SaveNeighbourData();	//Again Save Nb data
	
	cout << "SM Pairs" <<endl;
	for (int k = 0 ; k< dom.Nproc; k++)
		cout << dom.SMPairs[k].Size() <<endl;
	
	// cout << "Contact Pairs" <<endl;
	// for (int k = 0 ; k< dom.Nproc; k++)
		// cout << dom.ContPairs[k].Size() <<endl;

	cout << "FSM Pairs after contact" <<endl;
	for (int k = 0 ; k< dom.Nproc; k++)
		cout << dom.FSMPairs[k].Size() <<endl;
	
	dom.WriteXDMF("ContactTest");
}