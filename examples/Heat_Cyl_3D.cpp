
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
#include "InteractionAlt.cpp"
#include "SolverKickDrift.cpp"
#include "SolverVerlet.cpp"
#include "SolverFraser.cpp"

#define TAU		0.005
#define VMAX	10.0

int forcepart_count;

std::ofstream of;
double tout;

void UserAcc(SPH::Domain & domi)
{

}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
  SPH::Domain	dom;

  dom.Dimension	= 3;
  dom.Nproc	= 8;
  dom.Kernel_Set(Qubic_Spline);
  //dom.Kernel_Set(Hyperbolic_Spline);
  dom.Scheme	= 1;	//Mod Verlet
  //dom.XSPH	= 0.1; //Very important

  double dx,h,rho,K,G,Cs,Fy;
  double R,L,n;

  R	= 0.09;
  L	= 0.35;
  n	= 30.0;		//in length, radius is same distance

  rho	= 8410.0;
  K	= 6.7549e10;
  G	= 2.5902e10;
  Fy	= 300.e6;
  //dx	= L / (n-1);
  //dx = L/(n-1);
  dx = 0.008;
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
  
  //dom.Gradient_Approach_Set( 2); 


		// inline void Domain::AddCylinderLength(int tag, Vec3_t const & V, double Rxy, double Lz, 
									// double r, double Density, double h, bool Fixed) {
										
		//dom.AddCylinderLength(0, Vec3_t(0.,0.,0.0), R, L ,  dx/2., rho, h, false); 
		
    // void Domain::AddCylUniformLength(int tag, double Rxy, double Lz, 
																				// double r, double Density, double h, double ang, int rows, double r_i, bool ghost) {
                                          
    //dom.AddCylinderLength(0, Vec3_t(0.,0.,0.0), R, L ,  dx/2., rho, h, false); 
    
    dom.AddCylUniformLength(0, R, L ,  dx/2., rho, h, M_PI ,0, 0.0, false); 

    
		cout << "Particle count: "<<dom.Particles.Size()<<endl;
		
		forcepart_count = 0;
    int bottom_count = 0;
		//dom.gradKernelCorr = true;
		dom.ts_nb_inc = 5;		
		
      double x, y;
    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		x = dom.Particles[a]->x(0);
    		y = dom.Particles[a]->x(1);
        
			dom.Particles[a]->k_T			=	160.0;
			dom.Particles[a]->cp_T			=	560.0;
			dom.Particles[a]->T				= 700.0;		
      
    		dom.Particles[a]->G		= G;
    		dom.Particles[a]->PresEq	= 0;
    		dom.Particles[a]->Cs		= Cs;
    		dom.Particles[a]->Shepard	= false;
    		dom.Particles[a]->Material	= 2;
    		dom.Particles[a]->Fail		= 1;
          dom.Particles[a]->T_inf = 22.0;
          
    		//dom.Particles[a]->Beta		= 1.0;
    		dom.Particles[a]->TI		= 0.3;
    		dom.Particles[a]->TIInitDist	= dx;
    		double z = dom.Particles[a]->x(2);
    		//if ( z < -3.0*dx ){
        if ( sqrt(x*x+y*y)  > (R-dx) || z > (L-dx) || z < dx){
    			dom.Particles[a]->ID=2;
    			dom.Particles[a]->Thermal_BC 	= TH_BC_CONVECTION;
          dom.Particles[a]->T_inf = 22.0;
          dom.Particles[a]->h_conv		= 500.0; //W/m2-K
	    			// dom.Particles[a]->IsFree=false;
    			// dom.Particles[a]->NoSlip=true;			
          bottom_count++;
				}

    	}

      cout << "Conv  count "<<bottom_count<<endl;
		dom.WriteXDMF("maz");
		dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
		dom.BC.InOutFlow = 0;

    of = std::ofstream ("cf.csv", std::ios::out);
    of << "Time, disp, cf, ext_f_wk, plastic_wk, heat_cond, friction_wk"<<endl;
    tout = 0.;
    
    //dom.Solve_orig_Ext(/*tf*/0.00205,/*dt*/timestep,/*dtOut*/0.001,"test06",999);
		//dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/0.0001,"test06",999);
    timestep = (0.3*h*h*rho*dom.Particles[0]->cp_T/dom.Particles[0]->k_T);	
    timestep = 0.001;
    dom.thermal_solver = true;
        
    dom.CFL = 1.0;
    //timestep = 2.5e-6;
    dom.auto_ts = false;
    //dom.SolveDiffUpdateKickDrift(/*tf*/0.105,/*dt*/timestep,/*dtOut*/1.e-4,"test06",10000);	
    //dom.SolveDiffUpdateLeapfrog(/*tf*/0.105,/*dt*/timestep,/*dtOut*/1.e-4,"test06",10000);	
    //dom.SolveDiffUpdateVerlet(/*tf*/0.105,/*dt*/timestep,/*dtOut*/1.e-4,"test06",10000);	
    dom.ThermalSolve(/*tf*/10.0,/*dt*/timestep,/*dtOut*/1.e-1,"test06",10000);	
    //dom.SolveDiffUpdateFraser(5*timestep,/*dt*/timestep,/*dtOut*/timestep,"test06",10000);	
  
  
		return 0;
}
MECHSYS_CATCH
