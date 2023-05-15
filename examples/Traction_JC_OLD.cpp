// // **********************************************************************************
// // * PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids        * 
// // *             and soils) using Smoothed Particle Hydrodynamics method              *   
// // * Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo-Torres              *
// // *                                                                                  *
// // * This file is part of PersianSPH                                                  *
// // *                                                                                  *
// // * This is free software; you can redistribute it and/or modify it under the        *
// // * terms of the GNU General Public License as published by the Free Software        *
// // * Foundation; either version 3 of the License, or (at your option) any later       *
// // * version.                                                                         *
// // *                                                                                  *
// // * This program is distributed in the hope that it will be useful, but WITHOUT ANY  *
// // * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  *
// // * PARTICULAR PURPOSE. See the GNU General Public License for more details.         *
// // *                                                                                  *
// // * You should have received a copy of the GNU General Public License along with     *
// // * PersianSPH; if not, see <http://www.gnu.org/licenses/>                           *
// // ************************************************************************************/

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
  
  //vtraction = 0.1; //DIviding by Lo of Traction probe gives str_rate=1.
	
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

int main(int argc, char **argv) try
{
      SPH::Domain	dom;

      dom.Dimension	= 3;
      dom.Nproc	= 4;
    	dom.Kernel_Set(Qubic_Spline);
    	dom.Scheme	= 1;	//Mod Verlet
			//dom.XSPH	= 0.5; //Very important
                      

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

      Elastic_ el(E,nu);              

///// MATERIAL CONSTANTS EXAMPLE FROM
///// Zhang_2017 (Aluminium)
      double A,B,C,n_,m,T_m,T_t,eps_0;
      A = 175.e6; B = 380.0e6; C = 0.0015;
      m = 1.0;  n_ = 0.34; eps_0 = 1.0;
      T_m = 775.; T_t = 273.;
      
//Material A (MPa) B (MPa) C      n     m   Troom (K) Tmelt (K)
//Aluminum 175      380    0.0015 0.34 1.00 273       775			
//JohnsonCook(const double &a, const double &b, const double &c, const double &eps_0):
	// JohnsonCook(const double &a, const double &b, const double &c, 
              // const double &m_, const double &n_, const double &eps_0_, 
              // const double &T_m_, const double &T_t_)
			// JohnsonCook mat(el, A,B,C,
                      // m,n_,eps_0,
                      // T_m, T_t);
											
		Hollomon mat(el,Fy/E,1220.e6,0.195);
                      
    	rho	= 7850.0;
		K= E / ( 3.*(1.-2*nu) );
		G= E / (2.* (1.+nu));
		Fy	= 350.e6;


		dx = 0.010;
    h	= dx*1.2; //Very important

        Cs	= sqrt(K/rho);

        double timestep;
        //timestep = (0.2*h/(Cs));
		
		//timestep = 2.5e-6;
		timestep = 5.e-7;

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
		//dom.auto_ts = false;
		dom.AddTractionProbeLength(1, Vec3_t(0.,0.,-Lz_side/5.), R, Lz_side + Lz_side/5.,
											Lz_neckmin,Lz_necktot,Rxy_center,
											dx/2., rho, h, false);


		cout << "Particle count: "<<dom.Particles.Size()<<endl;
		
		//Center particle id
		//4081 if dx=10mm
		//2421 if dx=12mm
		int id=4081; 
		//dom.Particles[6777]->print_history = true;
		dom.Particles[id]->print_history = true; //If dx=10mm
		dom.Particles[id]->ID = 4;
		//dom.Particles[2421]->print_history = true;	//Particle 2421, 3 [     -0.006    -0.006     0.242 ]		
		cout << "Initial Yield Stress: "<<Fy<<endl;
    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
				
				dom.Particles[a]-> Material_model = HOLLOMON;
			 dom.Particles[a]->mat = &mat;
        
        dom.Particles[a]-> T = 273.;
				// dom.Particles[a]->k_T			=	150.;
				// dom.Particles[a]->cp_T			=	960.; 
				
				dom.Particles[a]->G		= G;
    		dom.Particles[a]->PresEq	= 0;
    		dom.Particles[a]->Cs		= Cs;
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
		dom.WriteXDMF("maz");
//		dom.m_kernel = SPH::iKernel(dom.Dimension,h);	


    	//dom.ThermalStructSolve(/*tf0.0105,/*dt*/timestep,/*dtOut*/0.0001,"test06",999);
      dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/0.0001,"test06",999);
        return 0;
}
MECHSYS_CATCH
