
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

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
  SPH::Domain	dom;

  dom.Dimension	= 3;
  dom.Nproc	= 4;
  dom.Kernel_Set(Qubic_Spline);
  dom.Scheme	= 1;	//Mod Verlet

	double dx,h,rho;
	double L;

	L	= 1.0;		
  rho	= 1.0;
  dx = 0.25;
  h	= dx*1.2; //Very important

  dom.DomMax(0) = L;
  dom.DomMin(0) = -L;
  
  Vec3_t p0 = Vec3_t ( 2.0, 0., 0.);
										
	//dom.AddCylinderLength(1, Vec3_t(0.,0.,-L/10.), R, L + 2.*L/10.,  dx/2., rho, h, false); 

// inline void Domain::AddBoxLength(int tag, Vec3_t const & V, double Lx, double Ly, double Lz, 
									// double r, double Density, double h, int type, int rotation, bool random, bool Fixed)
                  
  dom.AddBoxLength(1 ,p0, 
                      L + dx/10.0 , dx ,  dx , 
                      dx/2.0 ,rho, h, 1 , 0 , false, false );
      
  cout << "Particle count: "<<dom.Particles.Size()<<endl;
  dom.WriteXDMF("maz");
  
  //Write Kernels
  //A CORRECTIVE SMOOTHED PARTICLE METHOD FOR
  //BOUNDARY VALUE PROBLEMS IN HEAT CONDUCTION
  //Chen 1999
	for (int k=0; k<dom.Nproc;k++) {
    double mi,mj;
    Particle *P1,*P2;
    Vec3_t xij;
		for (size_t a=0; a<SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			//cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
			P1	= Particles[SMPairs[k][a].first];
			P2	= Particles[SMPairs[k][a].second];
			xij	= P1->x - P2->x;
						
			mi = P1->Mass;
			mj = P2->Mass;	
			//Eqn 3-112 Fraser Thesis
			P1->normal += mj * xij; 
			P2->normal -= mi * xij;					
		} //Nproc //Pairs  
  }
  //dom.m_kernel = SPH::iKernel(dom.Dimension,h);	

  return 0;
}
MECHSYS_CATCH
