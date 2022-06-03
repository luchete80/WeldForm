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
void UserAcc(SPH::Domain & domi)
{
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

/////////////// COMPACT NSEARCH ////////////////////
#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <chrono>
#include <random>

using namespace CompactNSearch;

void Plate_Al_Example(SPH::Domain &dom);

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

	ofstream outmesh; // outdata is like cin
	outmesh.open("outmesh.txt"); // opens the file
	
	cout << "Generating domain"<<endl;
	
	//THIS IS FOR A SIMPLE TEST
    dom.AddBoxLength(1 ,Vec3_t ( -L/2.0-L/20.0 , -H/2.0 , 0.0 ), L + L/10.0 + dx/10.0 , H + dx/10.0 ,  0 , dx/2.0 ,rho, h, 1 , 0 , false, false );
    
	for (int p=0;p < dom.Particles.Size();p++){
		 
		 outmesh << p << ", "<<dom.Particles[p]->x[0]<<", "<<dom.Particles[p]->x[1]<< ", " << dom.Particles[p]->x[2] <<endl;
	 }
	 outmesh.close();
	Real const radius =  static_cast<Real>(2.*h) ;	
    
	
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
	
	//
	// THIS IS NOT WORKING
	// //Pass it to SMPairs
	// //MAKING IT MULTITHREAD
	// // DYNAMICALLY ASSINGED!
	// it = neigbours_set.begin();
	// #pragma omp parallel for schedule (dynamic) num_threads(dom.Nproc)
	// for (int i=0;i<neigbours_set.size();i++){
		// size_t T = omp_get_thread_num();
				
			// dom.SMPairs[T].Push(std::make_pair(it->first, it->second));
		// it++;
	// }

	
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
	
	Plate_Al_Example(dom);

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

void Plate_Al_Example(SPH::Domain &dom){
	
        double dx,h,rho,K,G,Cs,Fy;
    	double H,L,n;	

        dom.Dimension	= 2;
        dom.Nproc	= 4;
    	dom.Kernel_Set(Quintic_Spline);
    	dom.Scheme	= 0;
		
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

    	for (size_t a=0; a<dom.Particles.Size(); a++)
    		dom.Particles[a]->h = h;	
			
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
    	dom.Solve_wo_init(/*tf*/0.01,/*dt*/timestep,/*dtOut*/0.001,"test06",999);
	
}