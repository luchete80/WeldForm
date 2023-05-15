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
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	for (int i=0; i<domi.Particles.Size(); i++)
	{
		if (domi.Particles[i]->ID == 3)
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(1.0e-2,0.0,0.0);
			domi.Particles[i]->vb		= Vec3_t(1.0e-2,0.0,0.0);
//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		}
		if (domi.Particles[i]->ID == 2)
		{
			domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			domi.Particles[i]->v		= Vec3_t(-1.0e-2,0.0,0.0);
			domi.Particles[i]->vb		= Vec3_t(-1.0e-2,0.0,0.0);
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
    	dom.Kernel_Set(Quintic_Spline);
    	dom.Scheme	= 0;
//     	dom.XSPH	= 0.5; //Very important

        double dx,h,rho,K,G,Cs,Fy;
    	double H,L,n;

    	H	= 0.01;
    	L	= 0.03;
    	n	= 15.0;	//ORIGINAL IS 40
		
		//lambda = E nu / ((1+nu)*(1-2 nu)) 
		//G = E/(2*(1+nu) )
		double E=72.e9;
		double nu=0.3;
		
    	rho	= 2800.0;
    	//K	= 3.25e6;
    	//G	= 7.15e5;
		K= E / ( 3.*(1.-2*nu) );
		G= E / (2.* (1.+nu));
		Fy	= 500.0e6;
    	dx	= H / n;
    	h	= dx*0.5; //Very important	//COMPARE WITH ANOTHER VALUES
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
		
     	dom.AddBoxLength(1 ,Vec3_t ( -L/2.0 -L/20.0, -H/2.0 , -H/2.0 ), 
							L + L/10.0 +dx/10.0 , H + dx/10.0 ,  H + dx/10.0 , 
							dx/2.0 ,rho, h, 1 , 0 , false, false );
		
		cout << "Particle count: "<<dom.Particles.Size()<<endl;
     	double x;

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		dom.Particles[a]->G			= G;
    		dom.Particles[a]->PresEq	= 0;
    		dom.Particles[a]->Cs		= Cs;
    		dom.Particles[a]->Shepard	= false;
    		dom.Particles[a]->Material	= 2;
    		dom.Particles[a]->Fail		= 1;
    		dom.Particles[a]->Sigmay	= Fy;
    		dom.Particles[a]->Alpha		= 1.0;
    		dom.Particles[a]->TI		= 0.3;
    		dom.Particles[a]->TIInitDist	= dx;
    		x = dom.Particles[a]->x(0);
    		if (x<-L/2.0)
    			dom.Particles[a]->ID=2;
    		if (x>L/2.0)
    			dom.Particles[a]->ID=3;
    	}

	
//    	dom.WriteXDMF("maz");
    	dom.Solve(/*tf*/0.05,/*dt*/timestep,/*dtOut*/1.e-5,"test06",100000);
        return 0;
}
MECHSYS_CATCH
