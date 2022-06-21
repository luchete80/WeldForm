#pragma once
#include <fstream>
#include <iostream>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include <iomanip>	//ONY FOR GCC!!


template <typename T>
bool readValue(const nlohmann::json &j, T &v)
{
	if (j.is_null())
		return false;

	v = j.get<T>();
	return true;
}


#include "Domain.h"

#define TAU		0.005
#define VMAX	1.0



void UserAcc(SPH::Domain & domi)
{
	double vtraction;

	if (domi.getTime() < TAU ) 
		vtraction = VMAX/TAU * domi.getTime();
	else
		vtraction = VMAX;
	
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)

	#ifdef __GNUC__
	for (size_t i=0; i<domi.Particles.Size(); i++)
	#else
	for (int i=0; i<domi.Particles.Size(); i++)
	#endif
	
	{
		if (domi.Particles[i]->ID == 3)
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,vtraction);
			domi.Particles[i]->va		= Vec3_t(0.0,0.0,vtraction);
			domi.Particles[i]->vb		= Vec3_t(0.0,0.0,vtraction);
//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
		if (domi.Particles[i]->ID == 2)
		{
			// domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			// domi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);
			// domi.Particles[i]->vb		= Vec3_t(0.0,0.0,0.0);
//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
	}
}


using std::cout;
using std::endl;


// int main(){

	// std::ifstream i("input.json");
	// json j;
	// i >> j;
	
	// //std::cout << j{"Configuration"}["pause"] << std::endl;    // returns boolean
	
	// nlohmann::json config = j["Configuration"];
	
	// double ts;
	// //config["timeStepSize"].get<ts>;
	// readValue(config["timeStepSize"], /*scene.timeStepSize*/ts);
	// std::cout << "Time Step size is: "<<ts<<std::endl;
	// // for(auto &array : j["objList"]) {
    // // std::cout << array["key1"] << std::endl;    // returns string
    // // std::cout << array["key2"] << std::endl;    // returns array
	// // }

	// // write prettified JSON to another file
	// std::ofstream o("pretty.json");
	// o << std::setw(4) << j << std::endl;
	

	
// }


int main(int argc, char **argv) try
{
	
	std::ifstream i("input.json");
	json j;
	i >> j;
	nlohmann::json config = j["Configuration"];	
  
	SPH::Domain	dom;

	dom.Dimension	= 3;
	dom.Nproc	= 4;
	dom.Kernel_Set(Qubic_Spline);
	dom.Scheme	= 1;	//Mod Verlet
	//dom.XSPH	= 0.5; //Very important

	double dx,K,G,Cs;
	double R,L,n;
	double Lz_side,Lz_neckmin,Lz_necktot,Rxy_center;
		
	R	= 0.075;

	Lz_side =0.2;
	Lz_neckmin = 0.050;
	Lz_necktot = 0.100;
	Rxy_center = 0.050;
	L = 2. * Lz_side + Lz_necktot;

	double ts;
	//config["timeStepSize"].get<ts>;
	readValue(config["timeStepSize"], /*scene.timeStepSize*/ts);
	std::cout << "Time Step size is: "<<ts<<std::endl;

	double radius;
	readValue(config["particleRadius"], radius);
	dx = 2.0*radius;//0.0085;
	double hFactor,h;
	
	readValue(config["hFactor"], h);
	h	= dx*hFactor; //Very important
	
	config = j["Material"];	
		
	double E;//  = 210.e9;
	readValue(config["youngsModulus"], E);
	double nu;// = 0.3;
	readValue(config["poissonsRatio"], nu);
  double rho;//	= 7850.0;	
	readValue(config["density0"], rho);
	std::cout << "Young's Modulus is : "<<E<<std::endl;
		

	K= E / ( 3.*(1.-2*nu) );
	G= E / (2.* (1.+nu));

	double Fy;//	= 350.e6;
	readValue(config["yieldStress0"], Fy);
	
	Cs	= sqrt(K/rho);

	double timestep;
	timestep = (0.2*h/(Cs));
		
	cout<<"t  = "<<timestep<<endl;
	cout<<"Cs = "<<Cs<<endl;
	cout<<"K  = "<<K<<endl;
	cout<<"G  = "<<G<<endl;

	dom.GeneralAfter = & UserAcc;
	dom.DomMax(0) = L;
	dom.DomMin(0) = -L;


	// inline void Domain::AddCylinderLength(int tag, Vec3_t const & V, double Rxy, double Lz, 
						// double r, double Density, double h, bool Fixed) {

	dom.AddTractionProbeLength(1, Vec3_t(0.,0.,-Lz_side/5.), R, Lz_side + Lz_side/5.,
								Lz_neckmin,Lz_necktot,Rxy_center,
								dx/2., rho, h, false);


	cout << "Particle count: "<<dom.Particles.Size()<<endl;

	for (size_t a=0; a<dom.Particles.Size(); a++) {
		dom.Particles[a]->G				= G;
		dom.Particles[a]->PresEq	= 0;
		dom.Particles[a]->Cs			= Cs;
		dom.Particles[a]->Shepard	= false;
		dom.Particles[a]->Material	= 2;
		dom.Particles[a]->Fail		= 1;
		dom.Particles[a]->Sigmay	= Fy;
		dom.Particles[a]->Alpha		= 1.0; 
		dom.Particles[a]->TI			= 0.3;
		dom.Particles[a]->TIInitDist	= dx;
		double z = dom.Particles[a]->x(2);
		if ( z < 0 ){
			dom.Particles[a]->ID=2;
			dom.Particles[a]->IsFree=false;
			dom.Particles[a]->NoSlip=true;    		
		}
		if ( z > L )
			dom.Particles[a]->ID=3;
	}
	//dom.WriteXDMF("maz");


	dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/0.0001,"test06",999);
  return 0;
}
MECHSYS_CATCH
