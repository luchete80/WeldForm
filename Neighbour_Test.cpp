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

#define PRINT(x)	cout << x[0]<<", "<<x[1]<<", "<<x[2]<<endl;
void UserAcc(SPH::Domain & domi) {

}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
   SPH::Domain	dom;

	dom.Dimension	= 3;
	dom.Nproc	= 4;
	dom.Kernel_Set(Quintic_Spline);
	dom.Scheme	= 0;
//     	dom.XSPH	= 0.5; //Very important
	double H,L,n,dx;

	H	= 0.01;
	L	= 0.03;
	n	= 30.0;	//ORIGINAL IS 40
	
	dx	= H / n;
	double h	= dx*1.3; //Very important

	dom.GeneralAfter = & UserAcc;
	dom.DomMax(0) = H;
	dom.DomMin(0) = -H;
	double rho	= 1000.0;

	dom.AddBoxLength(1 ,Vec3_t ( -H/2.0 , -H/2.0 , 0.0 ), H , H ,  H  , dx/2.0 ,rho, h, 1 , 0 , false, false );
	
	cout << "Particle count: "<<dom.Particles.Size()<<endl;
	double x;

	//dom.InitialChecks();
	dom.CellInitiate();
	dom.ListGenerate();
	// PrintInput(TheFileKey);
	// TimestepCheck();
	// WholeVelocity();
	
	double clock_time_spent,acc_time_spent;
	double neigbour_time_spent_per_interval;
	clock_t clock_beg;
	unsigned long steps=0;
	unsigned int first_step;
	for (size_t a=0; a< 20 ; a++) {
		clock_beg=clock();
		dom.MainNeighbourSearch();
		double t=(double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		neigbour_time_spent_per_interval += t;
		cout << "Total Neighbour search time in this interval: " <<  t <<endl;
		// for (int p=0;p<dom.Particles.size();p++){
			// cout << p <<": ";PRINT(dom.Particles[p]->x);
		// }

		dom.CellReset();
		dom.ListGenerate();
		steps++;
	}
		cout << "Average Neighbour search time in this interval: " << neigbour_time_spent_per_interval/(float)(steps)<<endl;

    	dom.WriteXDMF("maz");

	return 0;
}
MECHSYS_CATCH
