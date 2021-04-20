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

#include <iostream>
using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;


#define PRINT(x)	cout << x[0]<<", "<<x[1]<<", "<<x[2]<<endl;
void UserAcc(SPH::Domain & domi) {

}

/////////////// COMPACT NSEARCH ////////////////////
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <chrono>
#include <random>

using namespace CompactNSearch;



std::size_t const N_enright_steps = 50;

Real
compute_average_number_of_neighbors(NeighborhoodSearch const& nsearch);

Real
compute_average_distance(NeighborhoodSearch const& nsearch);

/////////////////////////////

#include <vector>;

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
  SPH::Domain	dom;

	dom.Dimension	= 3;
	dom.Nproc	= 1;
	dom.Kernel_Set(Quintic_Spline);
	dom.Scheme	= 1;
//     	dom.XSPH	= 0.5; //Very important
	double H,L,n,dx;

	H	= 0.01;
	L	= 0.03;
	n	= 15.0;	//ORIGINAL IS 40
	
	dx	= H / n;
	double h	= dx*0.5; //Very important

	dom.GeneralAfter = & UserAcc;
	dom.DomMax(0) = H;
	dom.DomMin(0) = -H;
	double rho	= 1000.0;

	ofstream outmesh; // outdata is like cin
	outmesh.open("outmesh.txt"); // opens the file
	
	cout << "Generating domain"<<endl;
	
	//dom.AddBoxLength(1 ,Vec3_t ( -H/2.0 , -H/2.0 , -H/2.0 ), H , H ,  H  , dx/2.0 ,rho, h, 1 , 0 , false, false );
     	dom.AddBoxLength(1 ,Vec3_t ( -L/2.0 , -H/2.0 , -H/2.0 ), 
							L + dx/10.0 , H + dx/10.0 ,  H + dx/10.0 , 
							dx/2.0 ,rho, h, 1 , 0 , false, false );
     for (int p=0;p < dom.Particles.Size();p++){
		 
		 outmesh << p << ", "<<dom.Particles[p]->x[0]<<", "<<dom.Particles[p]->x[1]<< ", " << dom.Particles[p]->x[2] <<endl;
	 }
	 outmesh.close();
	Real const r_omega = static_cast<Real>(H/2.)/ static_cast<Real>(n - 1);
	Real const radius =  static_cast<Real>(2.2) * static_cast<Real>(2.) * r_omega;	
    
	
	//dom.WriteXDMF("maz");
	///////////////////////////// COMPACT SEARCH THING
	std::vector<std::array<Real, 3>> positions;

	// for (unsigned int i = 0; i < (3*n); ++i){
		// for (unsigned int j = 0; j < n; ++j){
			// for (unsigned int k = 0; k < n; ++k) {
				// // std::array<Real, 3> x = {{
						
						// // r_omega * (static_cast<Real>(2.0) * static_cast<Real>(i) ),
						// // r_omega * (static_cast<Real>(2.0) * static_cast<Real>(j) ),
						// // r_omega * (static_cast<Real>(2.0) * static_cast<Real>(k) )						
						
						// // }};
				
				// std::array<Real, 3> x = {{
					// static_cast<Real>(dom.Particles[p]->x[0]),
					// static_cast<Real>(dom.Particles[p]->x[1]),
					// static_cast<Real>(dom.Particles[p]->x[2])
				// }};
				// positions.push_back(x);

			// }
		// }
	// }

	for (unsigned int p = 0; p < dom.Particles.Size(); p++){
		
		std::array<Real, 3> x = {{
			static_cast<Real>(dom.Particles[p]->x[0]),
			static_cast<Real>(dom.Particles[p]->x[1]),
			static_cast<Real>(dom.Particles[p]->x[2])
		}};
		positions.push_back(x);

	}
	

	
	//std::random_shuffle(positions.begin(), positions.end());
	cout << "Finding NEighbours"<<endl;
	NeighborhoodSearch nsearch(radius, true);
	nsearch.add_point_set(positions.front().data(), positions.size(), true, true);
	nsearch.add_point_set(positions.front().data(), positions.size(), true, true);
	nsearch.find_neighbors();

	nsearch.update_point_sets();
	std::vector<std::vector<unsigned int>> neighbors2;
	nsearch.find_neighbors(0, 1, neighbors2);
	std::vector<std::vector<unsigned int>> neighbors3;
	nsearch.find_neighbors(1, 2, neighbors3);


	//Pass to domain
	//std::set< std:: pair<int,int> > neigbours_set;
	auto const& d = nsearch.point_set(0);
	for (int i = 0; i < d.n_points(); ++i){
		const std::vector<unsigned int>& nbs = d.neighbor_list(0, i);
		//res += static_cast<unsigned long>(d.n_neighbors(0, i));
		for (int k=0;k< d.n_neighbors(0, i);k++) {
			//int n = d.neighbor(,i,k);
			
		}
	}


	std::cout << "#Points                                = " << positions.size() << std::endl;
	std::cout << "Search radius                          = " << radius << std::endl;
	// std::cout << "Min x                                  = " << min_x << std::endl;
	// std::cout << "Max x                                  = " << max_x << std::endl;
	std::cout << "Average number of neighbors            = " << compute_average_number_of_neighbors(nsearch) << std::endl;
	std::cout << "Average index distance prior to z-sort = " << compute_average_distance(nsearch) << std::endl;

	nsearch.z_sort();
	for (auto i = 0u; i < nsearch.n_point_sets(); ++i)
	{
		auto const& d = nsearch.point_set(i);
		d.sort_field(positions.data());

	}
	nsearch.find_neighbors();

	//compare_single_query_with_bruteforce_search(nsearch);
	//compare_with_bruteforce_search(nsearch);

	std::cout << "Average index distance after z-sort    = " << compute_average_distance(nsearch) << std::endl;

	std::cout << "Moving points:" << std::endl;
	for (int i = 0; i < N_enright_steps; ++i)
	{
		std::cout << "Enright step " << i << ". ";
		//advect();
		auto t0 = std::chrono::high_resolution_clock::now();
		nsearch.find_neighbors();
		std::cout << "Neighborhood search took " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count() << "ms" << std::endl;
		//compare_with_bruteforce_search(nsearch);
		//compare_single_query_with_bruteforce_search(nsearch);
	}
	

	return 0;
}
MECHSYS_CATCH

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

Real
compute_average_distance(NeighborhoodSearch const& nsearch)
{
	unsigned long long res = 0;
	auto const& d = nsearch.point_set(0);
	unsigned long long count = 0;
	for (int i = 0; i < d.n_points(); ++i)
	{
		std::size_t nn = d.n_neighbors(0, i);
		for (int j = 0; j < nn; ++j)
		{
			unsigned int k = d.neighbor(0, i, j);
			res += std::abs(i - static_cast<int>(k));
			count++;
		}
	}
	return static_cast<Real>(res) / static_cast<Real>(count);
}