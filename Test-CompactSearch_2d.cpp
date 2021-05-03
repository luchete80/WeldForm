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

#include <set>

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

void Test_Neigh();

std::size_t const N_enright_steps = 10;

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

	dom.Dimension	= 2;
	dom.Nproc	= 4;
	dom.Kernel_Set(Quintic_Spline);
	dom.Scheme	= 1;
//     	dom.XSPH	= 0.5; //Very important
	double H,L,n,dx;

	H	= 0.01;
	L	= 0.03;
	n	= 15.0;	//ORIGINAL IS 40
	
	dx	= H / n;
	double h	= dx*1.1; //Very important

	dom.GeneralAfter = & UserAcc;
	dom.DomMax(0) = H;
	dom.DomMin(0) = -H;
	double rho	= 2800.0;

	ofstream outmesh; // outdata is like cin
	outmesh.open("outmesh.txt"); // opens the file
	
	cout << "Generating domain"<<endl;
	
	//THIS IS FOR A SIMPLE TEST
     	dom.AddBoxLength(1 ,Vec3_t ( -L/2.0-L/20.0 , -H/2.0 , 0.0 ), L + L/10.0 + dx/10.0 , H + dx/10.0 ,  0 , 
						dx/2.0 ,rho, h, 1 , 0 , false, false );
		for (int p=0;p < dom.Particles.Size();p++){
		 
		 outmesh << p << ", "<<dom.Particles[p]->x[0]<<", "<<dom.Particles[p]->x[1]<< ", " << dom.Particles[p]->x[2] <<endl;
	 }
	 outmesh.close();
	Real const r_omega = static_cast<Real>(H/2.)/ static_cast<Real>(n - 1);
	Real const radius =  static_cast<Real>(2.0) * static_cast<Real>(2.) * r_omega;	
    
	
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
	
	outfind2.close();
	outfind2.open("find2_set.txt"); // opens the file
	std::set<std:: pair<int,int>>::iterator it = neigbours_set.begin();
	for (int i=0;i<neigbours_set.size();i++){
		
	outfind2<<it->first<<", "<<it->second<<endl;
	it++;
	}	
	outfind2.close();

	dom.CellInitiate();
	cout << "SM Pairs Size" << dom.SMPairs.Size()<<endl;
	//Array<std::pair<size_t,size_t> >		Initial;
	
	// //Pass it to SMPairs
	// //MAKING IT MULTITHREAD
	// // DYNAMICALLY ASSINGED!
	it = neigbours_set.begin();
	#pragma omp parallel for schedule (dynamic) num_threads(dom.Nproc)
	for (int i=0;i<neigbours_set.size();i++){
		size_t T = omp_get_thread_num();
				
			dom.SMPairs[T].Push(std::make_pair(it->first, it->second));
		it++;
	}

	
	for (int p=0;p<dom.Nproc;p++){
		cout << "Processor "<< p << ", " << dom.SMPairs[p].size()<< " pairs" << endl;
		
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

	for (int i=0 ; i<dom.Nproc ; i++) 
		dom.SMPairs[i].Clear();
	
	Test_Neigh();

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

void Test_Neigh(){
  SPH::Domain	dom;

	dom.Dimension	= 2;
	dom.Nproc	= 4;
	dom.Kernel_Set(Quintic_Spline);
	dom.Scheme	= 1;
//     	dom.XSPH	= 0.5; //Very important
	double H,L,n,dx;

	H	= 0.01;
	L	= 0.03;
	n	= 15.0;	//ORIGINAL IS 40
	
	dx	= H / n;
	double h	= dx*0.35; //Very important

	dom.GeneralAfter = & UserAcc;
	dom.DomMax(0) = H;
	dom.DomMin(0) = -H;
	double rho	= 2800.0;

	ofstream outmesh; // outdata is like cin
	outmesh.open("outmesh.txt"); // opens the file
	
	cout << "Generating domain"<<endl;
	
	//THIS IS FOR A SIMPLE TEST
     	dom.AddBoxLength(1 ,Vec3_t ( -L/2.0-L/20.0 , -H/2.0 , 0.0 ), L + L/10.0 + dx/10.0 , H + dx/10.0 ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );
		
		
							
	dom.CellInitiate();
	dom.ListGenerate();
	
	cout << "Particles: "<<dom.Particles.Size()<<endl;

		
	for (int i = 0; i < N_enright_steps; ++i)
	{
		std::cout << "Enright step " << i << ". ";
		//advect();
		auto t0 = std::chrono::high_resolution_clock::now();
		dom.MainNeighbourSearch();
		if (i == 0 )
			for (int p=0;p<dom.Nproc;p++)
				cout << "Processor "<< p << ", " << dom.SMPairs[p].size()<< " pairs" << endl;	
			
	for (int i=0 ; i<dom.Nproc ; i++) 
		dom.SMPairs[i].Clear();
		dom.CellReset();
		dom.ListGenerate();	
	
		std::cout << "Neighborhood search took " << 
			std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count() << 
			"ms" << std::endl;

	}


	dom.MainNeighbourSearch();
	ofstream outfind2; // outdata is like cin
	outfind2.open("findorig_set.txt"); // opens the file	
	for (int p=0;p<dom.Nproc;p++){
		outfind2 << "Processor "<< p << ", " << dom.SMPairs[p].size()<< " pairs" << endl;
		for (int i=0;i<dom.SMPairs[p].size();i++){
			outfind2 << dom.SMPairs[p][i].first << ", " << dom.SMPairs[p][i].second << endl;
			
		}
	}
	
	std::vector <int> nb(dom.Particles.Size());
	std::vector <int> nbcount(dom.Particles.Size());
	for ( size_t k = 0; k < dom.Nproc ; k++) {
		for (size_t a=0; a<dom.SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
		//cout << "a: " << a << "p1: " << dom.SMPairs[k][a].first << ", p2: "<< dom.SMPairs[k][a].second<<endl;
			nb[dom.SMPairs[k][a].first ]+=1;
			nb[dom.SMPairs[k][a].second]+=1;			
		}
	}	
	unsigned int avg=0;
	for (int i=0;i<nb.size();i++){
		//cout << "Neigbour "<< i <<": "<<nb[i]<<endl;
		avg+=nb[i];
	}
	avg/=dom.Particles.Size();
	
	cout << "For h: "<< h << "Avg Neighbour count is: "<<avg<<endl;
	outfind2.close();
}