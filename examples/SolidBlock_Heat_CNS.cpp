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
#include <CompactNSearch>
#include "Domain.h"

void UserAcc(SPH::Domain & domi)
{
	// #pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	// for (size_t i=0; i<domi.Particles.Size(); i++)
	// {
		// if (domi.Particles[i]->ID == 3)
		// {
			// domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			// domi.Particles[i]->v		= Vec3_t(1.0e-2,0.0,0.0);
			// domi.Particles[i]->vb		= Vec3_t(1.0e-2,0.0,0.0);
// //			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		// }
		// if (domi.Particles[i]->ID == 2)
		// {
			// domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
			// domi.Particles[i]->v		= Vec3_t(-1.0e-2,0.0,0.0);
			// domi.Particles[i]->vb		= Vec3_t(-1.0e-2,0.0,0.0);
// //			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);
		// }
	// }
}

/////////////// COMPACT NSEARCH ////////////////////
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <chrono>
#include <random>
using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;

#include <set>

using namespace CompactNSearch;


using std::cout;
using std::endl;

Real
compute_average_number_of_neighbors(NeighborhoodSearch const& nsearch)
{
	unsigned long res = 0;
	auto const& d = nsearch.point_set(0);
	for (int i = 0; i < d.n_points(); ++i)
	{
		res += static_cast<unsigned long>(d.n_neighbors(0, i));
	}
	return static_cast<Real>(res) / static_cast<Real>(d.n_points());
}

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

    	H	= 1.;
    	n	= 15.0;

    	rho	= 1000.0;
    	dx	= H / n;
    	h	= dx*1.2; //Very important
        Cs	= sqrt(K/rho);

        double timestep;

        cout<<"t  = "<<timestep<<endl;
        cout<<"Cs = "<<Cs<<endl;
        cout<<"K  = "<<K<<endl;
        cout<<"G  = "<<G<<endl;
        cout<<"Fy = "<<Fy<<endl;
    	dom.GeneralAfter = & UserAcc;
        dom.DomMax(0) = H;
        dom.DomMin(0) = -H;

     	dom.AddBoxLength(1 ,Vec3_t ( -H/2.0 -H/20., -H/2.0 -H/20., -H/2.0 -H/20. ), H + H/20., H +H/20.,  H + H/20. , dx/2.0 ,rho, h, 1 , 0 , false, false );

		Real const radius =  static_cast<Real>(2.*h);	
		
		
		//dom.WriteXDMF("maz");
		///////////////////////////// COMPACT SEARCH THING
		std::vector<std::array<Real, 3>> positions;

		for (unsigned int p = 0; p < dom.Particles.Size(); p++){
			
			std::array<Real, 3> x = {{
				static_cast<Real>(dom.Particles[p]->x[0]),
				static_cast<Real>(dom.Particles[p]->x[1]),
				static_cast<Real>(dom.Particles[p]->x[2])
			}};
			positions.push_back(x);

		}
		
	NeighborhoodSearch nsearch(radius, true);
	nsearch.add_point_set(positions.front().data(), positions.size(), true, true);
	nsearch.add_point_set(positions.front().data(), positions.size(), true, true);
	nsearch.find_neighbors();

	nsearch.update_point_sets();
	std::vector<std::vector<unsigned int>> neighbors2;
	nsearch.find_neighbors(0, 1, neighbors2);
	std::vector<std::vector<unsigned int>> neighbors3;
	nsearch.find_neighbors(1, 2, neighbors3);

	std::cout << "#Points                                = " << positions.size() << std::endl;
	std::cout << "Search radius                          = " << radius << std::endl;
	// std::cout << "Min x                                  = " << min_x << std::endl;
	// std::cout << "Max x                                  = " << max_x << std::endl;
	std::cout << "Average number of neighbors            = " << compute_average_number_of_neighbors(nsearch) << std::endl;
	

	ofstream outfind2; // outdata is like cin
	outfind2.open("find2_part.txt"); // opens the file
	//Pass to domain
	std::set< std:: pair<int,int> > neigbours_set;
	auto const& d = nsearch.point_set(0);
	for (int i = 0; i < d.n_points(); ++i){
		const std::vector<unsigned int>& nbs = d.neighbor_list(0, i);
		//res += static_cast<unsigned long>(d.n_neighbors(0, i));
		

		for (int k=0;k< nbs.size();k++) {
			outfind2<< i << ", "<<nbs[k]<<endl;
			neigbours_set.insert(std::make_pair( std::min(i,int(nbs[k])), std::max(i,int(nbs[k]))) );
			
		}
	}	
	cout << "Pares: "<<neigbours_set.size()<<endl;
	outfind2.close();

	std::set<std:: pair<int,int>>::iterator it = neigbours_set.begin();
	
	// outfind2.open("find2_set.txt"); // opens the file
	// for (int i=0;i<neigbours_set.size();i++){
		
	// outfind2<<it->first<<", "<<it->second<<endl;
	// it++;
	// }	
	// outfind2.close();
	
	//////////////////////////////
	////////CELL INITIATE IS CRUCIAL FOR INITIALIZE ARAYS!!!
	dom.CellInitiate();
	//////////////////////////////
	it = neigbours_set.begin();
	int pairsperproc = neigbours_set.size()/dom.Nproc;
	cout << "Pairs per proc: " <<pairsperproc<<endl;
	int pair=0;
	int nproc=0;
	while (it != neigbours_set.end()) {
		if (pair > (nproc + 1 ) * pairsperproc){
			nproc++;
			cout<<"changing proc"<< nproc<<", pair "<<pair<<endl;
		}
					
		dom.SMPairs[nproc].Push(std::make_pair(it->first, it->second));
		it++;
		pair++;
		
	}

	std::cout << "#Points                                = " << positions.size() << std::endl;
	std::cout << "Search radius                          = " << radius << std::endl;
	// std::cout << "Min x                                  = " << min_x << std::endl;
	// std::cout << "Max x                                  = " << max_x << std::endl;
	std::cout << "Average number of neighbors            = " << compute_average_number_of_neighbors(nsearch) << std::endl;
	
	
		std::cout << "Particle Number: "<< dom.Particles.size() << endl;
     	double x;

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		x = dom.Particles[a]->x(0);
			dom.Particles[a]->k_T			=	3000.;
			dom.Particles[a]->cp_T			=	1.;
			dom.Particles[a]->h_conv		= 100.0; //W/m2-K
			dom.Particles[a]->T_inf 		= 500.;
			dom.Particles[a]->T				= 20.0;			
    		if ( x < -H/2.0 ) {
    			dom.Particles[a]->Thermal_BC = TH_BC_CONVECTION;
				// cout << "Particle " << a << "is convection BC" <<endl;
			}
    	}

        timestep = (0.3*h*h*rho*dom.Particles[0]->cp_T/dom.Particles[0]->k_T);	
		cout << "Time Step: "<<timestep<<endl;
		//timestep=1.e-6;
		//0.3 rho cp h^2/k
	
	// std::vector <int> nb(dom.Particles.Size());
	// std::vector <int> nbcount(dom.Particles.Size());
	// for ( size_t k = 0; k < dom.Nproc ; k++) {
		// for (size_t a=0; a<dom.SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
		// //cout << "a: " << a << "p1: " << dom.SMPairs[k][a].first << ", p2: "<< dom.SMPairs[k][a].second<<endl;
			// nb[dom.SMPairs[k][a].first ]+=1;
			// nb[dom.SMPairs[k][a].second]+=1;
			
		// }
	// }	
	// for (int i=0;i<nb.size();i++)
		// cout << "Neigbour "<< i <<": "<<nb[i]<<endl;
		
//    	dom.WriteXDMF("maz");
//    	dom.Solve(/*tf*/0.01,/*dt*/timestep,/*dtOut*/0.001,"test06",999);

		dom.ThermalSolve_wo_init(/*tf*/1.,/*dt*/timestep,/*dtOut*/0.1,"test06",999);

//		dom.ThermalSolve(/*tf*/10.,/*dt*/timestep,/*dtOut*/0.1,"test06",999);

        return 0;
}


MECHSYS_CATCH
