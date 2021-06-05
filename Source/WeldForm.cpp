#include "Input.h"
/***********************************************************************************
* PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids        * 
*             and soils) using Smoothed Particle Hydrodynamics method              *   
* Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo-Torres              *
*                                                                                  *
* This file is part of PersianSPH                                                  *
*                                                                                  *
* This is free software; you can redistribute it and/or modify it under the        *
* terms of the GNU General Public License as published by the Free Software        *
* Foundation; either version 3 of the License, or (at your option) any later       *
* version.                                                                         *
*                                                                                  *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY  *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.         *
*                                                                                  *
* You should have received a copy of the GNU General Public License along with     *
* PersianSPH; if not, see <http://www.gnu.org/licenses/>                           *
************************************************************************************/

#include "Domain.h"
#include "Input.h"

#define TAU		0.005
#define VMAX	10.0

#define PRINTVEC(v)	cout << v[0]<<", "<<v[1]<<", "<<v[2]<<endl;

void UserAcc(SPH::Domain & domi)
{
	double vcompress;

	if (domi.getTime() < TAU ) 
		vcompress = VMAX/TAU * domi.getTime();
	else
		vcompress = VMAX;
	
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
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,-vcompress);

			if (domi.Scheme == 1 )
				domi.Particles[i]->va		= Vec3_t(0.0,0.0,-vcompress);//VERLET
			else
				domi.Particles[i]->vb		= Vec3_t(0.0,0.0,-vcompress);//LEAPFROG
//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
		if (domi.Particles[i]->ID == 2)
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,0.);
			if (domi.Scheme == 1 )
				domi.Particles[i]->va		= Vec3_t(0.0,0.0,-vcompress);//VERLET
			else
				domi.Particles[i]->vb		= Vec3_t(0.0,0.0,-vcompress);//LEAPFROG
//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
	}
}

using std::cout;
using std::endl;

int main(int argc, char **argv) try {

	if (argc > 1) {
		string inputFileName=argv[1];	
		std::ifstream i(argv[1]);
		json j;
		i >> j;
		
		nlohmann::json config 		= j["Configuration"];
		nlohmann::json material 	= j["Material"];
		nlohmann::json domblock 	= j["DomainBlock"];
		
		SPH::Domain	dom;
		
		bool sim2D;
		double ts;
		
		readValue(config["sim2D"], sim2D);
		
		if (sim2D)
			dom.Dimension	= 2;
        else
			dom.Dimension	= 3;
		
		dom.Nproc	= 4;
		string kernel;
		readValue(config["timeStepSize"], /*scene.timeStepSize*/ts);
    	dom.Kernel_Set(Qubic_Spline);
		
		readValue(config["integrationMethod"], dom.Scheme); //0:Verlet, 1:LeapFrog, 2: Modified Verlet

     	//dom.XSPH	= 0.5; //Very important

        double dx,r,h;


		readValue(config["particleRadius"], r);
		double hfactor;
		readValue(config["hFactor"], hfactor);

    	dom.GeneralAfter = & UserAcc;
		
	
		//////////////
		// MATERIAL //
		//////////////
		double rho,E,nu,K,G,Cs,Fy;
    	readValue(material["density0"], 		rho);
    	readValue(material["youngsModulus"], 	E);
    	readValue(material["poissonsRatio"], 	nu);
    	readValue(material["yieldStress0"], 	Fy);
		
		K= E / ( 3.*(1.-2*nu) );
		G= E / (2.* (1.+nu));

		dx 	= 2.*r;
    	h	= dx*hfactor; //Very important
        Cs	= sqrt(K/rho);

        double timestep,cflFactor;
		int cflMethod;
		readValue(config["cflMethod"], cflMethod);
		if (cflMethod == 0)
			readValue(config["timeStepSize"], timestep);
		else {
			readValue(config["cflFactor"], cflFactor);
			timestep = (cflFactor*h/(Cs));
		}

		////////////
		// DOMAIN //
		////////////
		Vec3_t start,L;
		readVector(domblock["start"], 	start);
		readVector(domblock["dim"], 	L);
        for (int i=0;i<3;i++) {//TODO: Increment by Start Vector
			dom.DomMax(0) = L[i];
			dom.DomMin(0) = -L[i];
		}		

		
		// inline void Domain::AddCylinderLength(int tag, Vec3_t const & V, double Rxy, double Lz, 
									// double r, double Density, double h, bool Fixed) {
												
		//dom.AddCylinderLength(1, Vec3_t(0.,0.,-L/10.), R, L + 2.*L/10. + dx, r, rho, h, false); 
		
		//dom.AddCylinderLength(1, start, L[0], L[2], r, rho, h, false); 

		

        cout <<"t  			= "<<timestep<<endl;
        cout <<"Cs 			= "<<Cs<<endl;
        cout <<"K  			= "<<E<<endl;
        cout <<"G  			= "<<nu<<endl;
        cout <<"Fy 			= "<<Fy<<endl;
		cout <<"dx 			= "<<dx<<endl;
		cout <<"h  			= "<<h<<endl;
		cout <<"-------------------------"<<endl;
		cout <<	"Dim: "<<dom.Dimension<<endl;				
		cout << "Particle count: "<<dom.Particles.Size()<<endl;

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		dom.Particles[a]->G				= G;
    		dom.Particles[a]->PresEq		= 0;
    		dom.Particles[a]->Cs			= Cs;
    		dom.Particles[a]->Shepard		= false;
    		dom.Particles[a]->Material		= 2;
    		dom.Particles[a]->Fail			= 1;
    		dom.Particles[a]->Sigmay		= Fy;
    		dom.Particles[a]->Alpha			= 0.0;
			dom.Particles[a]->Beta			= 0.0;
    		dom.Particles[a]->TI			= 0.3;
    		dom.Particles[a]->TIInitDist	= dx;
			
    		double z = dom.Particles[a]->x(2);
    		if ( z < 0 ){
    			dom.Particles[a]->ID=2;
				dom.Particles[a]->IsFree=false;
				dom.Particles[a]->NoSlip=true;
			} else if ( z > L[2] ){
    			dom.Particles[a]->ID=3;
				// dom.Particles[a]->IsFree=false;
				// dom.Particles[a]->NoSlip=true;
			}
    	}
		// dom.WriteXDMF("maz");
		// dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
		// dom.BC.InOutFlow = 0;

    	//dom.Solve(/*tf*/0.00205,/*dt*/timestep,/*dtOut*/0.00005,"test06",999);
	}	//Argc > 0
	
    return 0;
}

MECHSYS_CATCH