
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

#define TAU		0.005
#define VMAX	10.0



void UserAcc(SPH::Domain & domi)
{
	double vcompress;

	// if (domi.getTime() < TAU ) 
		// vcompress = VMAX/TAU * domi.getTime();
	// else
		// vcompress = VMAX;
	vcompress = 1.;
	//cout << "time: "<< domi.getTime() << "V compress "<< vcompress <<endl;
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
			domi.Particles[i]->va		= Vec3_t(0.0,0.0,-vcompress);
			//domi.Particles[i]->vb		= Vec3_t(0.0,0.0,-vcompress);
//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
		if (domi.Particles[i]->ID == 2)
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->va		= Vec3_t(0.0,0.0,0.0);
			//domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
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
		dx = 0.030;
    h	= dx*1.2; //Very important
        Cs	= sqrt(K/rho);

        double timestep;
        timestep = (0.2*h/(Cs));
				//timestep = 1.0e-7;
		
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
										
		dom.AddCylinderLength(1, Vec3_t(0.,0.,-L/10.), R, L + 2.*L/10.,  dx/2., rho, h, false); 
		
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
    		dom.Particles[a]->Alpha		= 0.0;
    		//dom.Particles[a]->Beta		= 1.0;
    		dom.Particles[a]->TI		= 0.0;
//    		dom.Particles[a]->TI		= 0.3;
    		dom.Particles[a]->TIInitDist	= dx;
    		double z = dom.Particles[a]->x(2);
    		if ( z < 0 ){
    			dom.Particles[a]->ID=2;
	    		//dom.Particles[a]->IsFree=false;
    			//dom.Particles[a]->NoSlip=true;			
				
				}
    		if ( z > L )
    			dom.Particles[a]->ID=3;
    	}
		dom.WriteXDMF("maz");
		dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
		dom.BC.InOutFlow = 0;

		dom.auto_ts = true;
		
		//IF performing a nb test	std::vector <int> nb(dom.Particles.Size());
		dom.CellInitiate();
		dom.ListGenerate();
		dom.MainNeighbourSearch();
		
		std::vector <int> nb(dom.Particles.Size());
		std::vector <int> nbcount(dom.Particles.Size());
		#pragma omp parallel for schedule (static) num_threads(dom.Nproc)
		for ( int k = 0; k < dom.Nproc ; k++) {
			for (int a=0; a<dom.SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			//cout << "a: " << a << "p1: " << dom.SMPairs[k][a].first << ", p2: "<< dom.SMPairs[k][a].second<<endl;
				if (dom.SMPairs[k][a].first == 0 )
					cout << dom.SMPairs[k][a].second <<", "<<endl;
				else if (dom.SMPairs[k][a].second == 0 )
					cout << dom.SMPairs[k][a].first <<", "<<endl;
			}
		}	
		//dom.Solve(/*tf*//*0.0105*/50.*timestep+1e-8,/*dt*/timestep,/*dtOut*/timestep,"test_gpu",999);
		dom.Solve(0.000101,/*dt*/timestep,/*dtOut*/0.00001,"test06",999);
    
		return 0;
}
MECHSYS_CATCH
