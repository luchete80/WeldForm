
#include "Domain.h"

#define TAU		0.005
#define VMAX	1.0


void UserAcc(SPH::Domain & domi)
{

}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
	SPH::Domain	dom;

	dom.Dimension	= 3;
	dom.Nproc	= 4;
	dom.Kernel_Set(Qubic_Spline);

	double dx,h,rho,K,G,Cs,Fy;
	double R,L,n;
	double Lz_side,Lz_neckmin,Lz_necktot,Rxy_center;

	R	= 0.075;

	Lz_side =0.2;
	Lz_neckmin = 0.050;
	Lz_necktot = 0.100;
	Rxy_center = 0.050;
	L = 2. * Lz_side + Lz_necktot;

	double E  = 210.e9;
	double nu = 0.3;

	rho	= 7850.0;
	K= E / ( 3.*(1.-2*nu) );
	G= E / (2.* (1.+nu));
	Fy	= 350.e6;


	dx = 0.008;
	h	= dx*1.1; //Very important

	Cs	= sqrt(K/rho);

	double timestep;
	//timestep = (0.2*h/(Cs));

	//timestep = 2.5e-6;
	timestep = 5.e-7;

	dom.GeneralAfter = & UserAcc;
	dom.DomMax(0) = L;
	dom.DomMin(0) = -L;

	dom.AddTractionProbeLength(1, Vec3_t(0.,0.,-Lz_side/10.), R, Lz_side + Lz_side/10.,
								Lz_neckmin,Lz_necktot,Rxy_center,
								dx/2., rho, h, false);


	cout << "Particle count: "<<dom.Particles.Size()<<endl;
	
	dom.CellInitiate();
	dom.ListGenerate();
	cout << "Main Nb Search"<<endl;
	dom.MainNeighbourSearch_Ext();	//MUST DO CELL INITIATE FIRST
	dom.SaveNeighbourData();
	
	
	dom.CalculateSurface(2);

	dom.WriteXDMF("maz");

	return 0;
}
MECHSYS_CATCH
