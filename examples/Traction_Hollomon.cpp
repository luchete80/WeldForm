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
#include "SolverFraser.cpp"

#define TAU		0.005
#define VMAX	1.0

#define DX 0.008
double tout, dtout;

ofstream ofprop("traction.csv", std::ios::out);
double L, dx;

void UserAcc(SPH::Domain & domi)
{
	double vtraction;

	if (domi.getTime() < TAU ) 
		vtraction = VMAX/TAU * domi.getTime();
	else
		vtraction = VMAX;

  double dS = DX * DX;
  double normal_acc_sum=0.;
  double max_seq = 0;
  int part = 0;
  for (size_t i=0; i<domi.Particles.Size(); i++){
    if (domi.Particles[i]->x(2) > (L -dx) && domi.Particles[i]->x(2) < L+dx/2.) {
      domi.m_scalar_prop  += domi.Particles[i]->Sigma (2,2) * dS;
      normal_acc_sum      += domi.Particles[i]->a(2) * domi.Particles[i]->Mass;
      part++;
      
    }
    
    if (domi.Particles[i]->Sigma_eq > max_seq)
      max_seq = domi.Particles[i]->Sigma_eq;
  }
  
  dtout = 1.0e-4;
  if (domi.getTime()>tout){
    cout << "Normal integrated force " <<domi.m_scalar_prop<<endl;
    cout << "Normal acc sum " << normal_acc_sum<<endl;
    tout += dtout; 
    cout << "particles acc "<<part<<endl;
    ofprop << domi.max_disp[0]<<", "<<domi.max_disp[1]<<", " <<domi.max_disp[2]<<", " << normal_acc_sum << ", " << max_seq<< endl;
  }
	
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
      //ext_work += domi.Particles[i]->Sigma (2,2) * DX * DX * domi.Particles[i]->Displacement(2);
//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
		if (domi.Particles[i]->ID == 2)
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->vb		= Vec3_t(0.0,0.0,0.0);
//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
	}
}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
      SPH::Domain	dom;
      tout = 0.;
      dom.Dimension	= 3;
      dom.Nproc	= 4;
    	dom.Kernel_Set(Qubic_Spline);
    	dom.Scheme	= 1;	//Mod Verlet
			//dom.XSPH	= 0.1; //Very important
			

        double h,rho,K,G,Cs,Fy;
    	double R,n;
		double Lz_side,Lz_neckmin,Lz_necktot,Rxy_center;
		
	R	= 0.075;

	Lz_side =0.2;
	Lz_neckmin = 0.050;
	Lz_necktot = 0.100;
	Rxy_center = 0.050;
	L = 2. * Lz_side + Lz_necktot;

	double E  = 200.e9;
	double nu = 0.3;

	rho	= 7850.0;
	K= E / ( 3.*(1.-2*nu) );
	G= E / (2.* (1.+nu));
	
	
	Elastic_ el(E,nu);
	//Hollomon(const double eps0_, const double &k_, const double &m_):
  /////ORIGINALLY 
  // Fy=350e6;
  // Hollomon mat(el,Fy,1220.e6,0.195);
	//////NEW
  Fy	= 260.e6;
  Hollomon mat(el,Fy,7.1568e8,0.22);
			

		dx = DX;
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
		//6777 if dx=8.5
		//4081 if dx=10mm
		//2421 if dx=12mm
		
		//dom.Particles[6777]->print_history = true;
		//dom.Particles[4081]->print_history = true;
		//dom.Particles[6777]->ID = 1;

		//dom.Particles[6777]->print_history = true;
		// dom.Particles[4081]->print_history = true; //If dx=10mm
		// dom.Particles[4081]->ID = 4;
		//If dx=12mm
		dom.Particles[4081]->print_history = true;	//Particle 2421, 3 [     -0.006    -0.006     0.242 ]	
    cout << "Particle History position: "<<endl;
    cout << dom.Particles[4081]->x(0)<<", "<<dom.Particles[4081]->x(1)<<", "<<dom.Particles[4081]->x(2)<<endl;
		int part_2, part_3;
    part_2 = part_3 = 0;
    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
				
				dom.Particles[a]-> Material_model = HOLLOMON;
				dom.Particles[a]->mat = &mat;
				
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
    		if ( z < dx ){
    			dom.Particles[a]->ID=2;
          part_2++;
    			// dom.Particles[a]->IsFree=false;
    			// dom.Particles[a]->NoSlip=true;    		
				}
				if ( z > L -dx){
    			dom.Particles[a]->ID=3;
          part_3++;
        }
    	}
		dom.WriteXDMF("maz");
    cout << "particles boundary"<<part_2<< ", "<<part_3<<endl;
//		dom.m_kernel = SPH::iKernel(dom.Dimension,h);	

  // // // dom.auto_ts=true;
  // // // dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/0.0001,"test06",999);

  timestep = (1.0*h/(Cs+VMAX)); //Standard modified Verlet do not accept such step
  //dom.auto_ts=false;
  dom.auto_ts=true;
  //dom.SolveDiffUpdateKickDrift(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/1.e-4 ,"test06",1000);
  
  ofprop << "maxux, maxux, maxux, force, maxseq "<<endl;
	dom.SolveDiffUpdateFraser(/*tf*/0.01005,/*dt*/timestep,/*dtOut*/1.e-4,"test06",1000);  
  return 0;
}
MECHSYS_CATCH
