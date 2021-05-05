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

void UserAcc(SPH::Domain & domi)
{

}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
       SPH::Domain	dom;

        dom.Dimension	= 2;
        dom.Nproc	= 4;
    	dom.Kernel_Set(Quintic_Spline);
    	dom.Scheme	= 0;
//     	dom.XSPH	= 0.5; //Very important

        double dx,h,rho,K,G,Cs,Fy;
    	double H,L,n;

    	H	= 0.01;
    	L	= 0.03;
    	n	= 40.0;	//ORIGINAL IS 40
		
    	rho	= 1000.0;
    	K	= 3.25e6;
    	G	= 7.15e5;
		Fy	= 4000.0;
    	dx	= H / n;
    	h	= dx*1.3; //Very important
        Cs	= sqrt(K/rho);

        double timestep;
        timestep = (0.2*h/(Cs));

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
										
		dom.AddCylinderLength(1, Vec3_t(0.,0.,0.), 0.01, 0.2, 
									5.e-4, 2800, h, false); 
		
		cout << "Particle count: "<<dom.Particles.Size()<<endl;

		dom.WriteXDMF("maz");

        return 0;
}
MECHSYS_CATCH
