
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

#include <vector>

double MyKernel(size_t const & Dim, size_t const & KT, double const & q, double const & h)
{
  double C;

  switch (KT)
  {
    case 0:	// Qubic Spline
      //Dim == 2 ? C = 10.0/(7.0*h*h*M_PI) : C = 1.0/(h*h*h*M_PI);
      C = 2./(3.*h);//DIM 1

      if 		  (q<1.0)	return C*(1.0-(3.0/2.0)*q*q+(3.0/4.0)*q*q*q);
      else if (q<2.0)	return C*((1.0/4.0)*(2.0-q)*(2.0-q)*(2.0-q));
      else						return 0.0;
      break;
  }
}

double MyGradKernel(size_t const & Dim, size_t const & KT, double const & q, double const & h)
{
  double C;

  switch (KT)
  {
    case 0:	// Qubic Spline
      //Dim ==2 ? C = 10.0/(7.0*h*h*h*M_PI) : C = 1.0/(h*h*h*h*M_PI);
      C = 2./(3.*h*h);
      if 		(q==0.0)	return C/h    *(-3.0+(9.0/2.0)*q);
      else if (q<1.0)		return C/(q*h)*(-3.0*q+(9.0/4.0)*q*q);
      else if (q<2.0)		return C/(q*h)*((-3.0/4.0)*(2.0-q)*(2.0-q));
      else							return 0.0;
      break;
  }
}
  
inline void MyListGenerate (SPH::Domain *dom)
{
	int i, temp=0;
  cout << "Cell size"<<dom->CellSize(0)<<endl;
  cout << "dom->BLPF(0)"<<dom->BLPF(0)<<endl; 
		for (size_t a=0; a<dom->Particles.Size(); a++)
		{
			i= (int) (floor((dom->Particles[a]->x(0) - dom->BLPF(0)) / dom->CellSize(0)));

			if (i<0)
            {
                    if ((dom->BLPF(0) - dom->Particles[a]->x(0)) <= dom->hmax) i=0;
                            else std::cout<<"Leaving i<0"<<std::endl;
            }

			if (i>=dom->CellNo[0])
			{
					if ((dom->Particles[a]->x(0) - dom->TRPR(0)) <= dom->hmax) i=dom->CellNo[0]-1;
							else std::cout<<"Leaving i>=CellNo"<<std::endl;
			}
      cout << "cell "<<i<<", part"<<a<<endl;
			temp = dom->HOC[i][0][0];
			dom->HOC[i][0][0] = a;
			dom->Particles[a]->LL = temp;
			dom->Particles[a]->CC[0] = i;
			dom->Particles[a]->CC[1] = 0;
			dom->Particles[a]->CC[2] = 0;
			if (!dom->Particles[a]->IsFree) dom->FixedParticles.Push(a);
    }
}

int main(int argc, char **argv) try
{
  SPH::Domain	dom;
  
  int Dimension = 2;
  
  dom.Dimension	= Dimension; // In fact is one
  dom.Nproc	= 1;
  dom.Kernel_Set(Qubic_Spline);
  dom.Scheme	= 1;	//Mod Verlet

	double dx,h,rho;
	double L;

	int n =5;
	L	= 1.0;		
  rho	= 1.0;
  dx = L/n;
  h	= dx*1.; //Very important
	
  dom.DomMax(0) = 2.;
  dom.DomMin(0) = 3.;
  
  Vec3_t p0 = Vec3_t ( 2.0, 0., 0.);
										
	//dom.AddCylinderLength(1, Vec3_t(0.,0.,-L/10.), R, L + 2.*L/10.,  dx/2., rho, h, false); 

// inline void Domain::AddBoxLength(int tag, Vec3_t const & V, double Lx, double Ly, double Lz, 
									// double r, double Density, double h, int type, int rotation, bool random, bool Fixed)
                  
  dom.AddBoxLength(1 ,p0 - dx/10, 
                      L + dx/10.0 , dx ,  0 , 
                      dx/2.0 ,rho, h, 1 , 0 , false, false );
      
  //cout << "Particle count: "<<dom.Particles.Size()<<endl;
  dom.WriteXDMF("maz");

  for (int i = 0; i<dom.Particles.Size();i++) 
    cout << "i" << i<< ", x: "<<dom.Particles[i]->x(0)<<endl;
    
    
  // cout << "Cell Init"<<endl;
	// dom.CellInitiate();
  // cout << "SMPairs size: "<<dom.SMPairs[0].Size()<<endl;
  // cout << "Generating list"<<endl;
  Array<std::pair<size_t,size_t> >           Initial;
  for(size_t i=0 ; i<dom.Nproc ; i++) 
		dom.SMPairs.Push(Initial);
	// MyListGenerate(&dom);
  // cout << "Searching for Nbs"<<endl;
	// dom.MainNeighbourSearch();
  // cout << "Done"<<endl;
  
  cout << "Inserting pairs"<<endl;

  dom.SMPairs[0].Push(std::make_pair(0, 1));
	for (int i=2;i<n;i++) {
		dom.SMPairs[0].Push(std::make_pair(i-1, i));
		dom.SMPairs[0].Push(std::make_pair(i-2, i));
	}

	
	std::vector <int> nb(n);
	
	//Case of n=5
  // dom.SMPairs[0].Push(std::make_pair(0, 1));
  // dom.SMPairs[0].Push(std::make_pair(0, 2));
  // dom.SMPairs[0].Push(std::make_pair(1, 2));
  // dom.SMPairs[0].Push(std::make_pair(1, 3));
  // dom.SMPairs[0].Push(std::make_pair(2, 3));
  // dom.SMPairs[0].Push(std::make_pair(2, 4));  
  // dom.SMPairs[0].Push(std::make_pair(3, 4));
  
  cout << "Done"<<endl;
  
  std::vector<double>  fx(dom.Particles.Size());
  std::vector<double> dfx(dom.Particles.Size());

  std::vector<double>  gx(dom.Particles.Size());
  //Write Kernels
  //A CORRECTIVE SMOOTHED PARTICLE METHOD FOR
  //BOUNDARY VALUE PROBLEMS IN HEAT CONDUCTION
  //Chen 1999

  cout << "Calculating Kernel..."<<endl;
	for (int k=0; k<dom.Nproc;k++) {
    double mi,mj;
    double di,dj;
    SPH::Particle *P1,*P2;
    Vec3_t xij;
		for (size_t a=0; a<dom.SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			//cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
      int i = dom.SMPairs[k][a].first;
      int j = dom.SMPairs[k][a].second;
			P1	= dom.Particles[i];
			P2	= dom.Particles[j];
			xij	= P1->x - P2->x;
						
			mi = P1->Mass;
			mj = P2->Mass;	
      
      di = P1->Density;
			dj = P2->Density;	
      double K	= MyKernel(Dimension, 0, norm(xij)/h, h);
      double GK = MyGradKernel(Dimension, 0, norm(xij)/h, h);
      
      cout <<"Vi"<<mj/dj<<endl;
      cout << "r, K, GK"<<norm(xij)/h<<", "<<K<<", "<<GK<<endl;
      fx[i] += /*mj/dj */dx * P2->x(0)*P2->x(0) * K;
      fx[j] += /*mi/di */dx * P1->x(0)*P1->x(0) * K;
      
      gx[i] += /*mj/dj*/ /*dx * */P2->x(0)* K;
      gx[j] += /*mi/di*/ /*dx * */P1->x(0)* K;
      
      dfx[i] -= dx * P2->x(0)*P2->x(0)*P2->x(0)/3 * GK * xij(0);
      dfx[j] += dx * P1->x(0)*P1->x(0)*P1->x(0)/3 * GK * xij(0);
			
			nb[i]++;
			nb[j]++;      
		} //Nproc //Pairs  
  }
	for (int i = 0; i<dom.Particles.Size();i++) {
		double x = dom.Particles[i]->x(0);
		double K	= MyKernel(Dimension, 0, 0, h);
		
		fx[i] += /*mj/dj */dx * x*x * K;

	}
  cout << "Done."<<endl;
  //dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
  
  for (int i = 0; i<dom.Particles.Size();i++) {
    cout << "Anal f" << dom.Particles[i]->x(0)<<", "<<dom.Particles[i]->x(0)*dom.Particles[i]->x(0)<<endl;
    cout << "x, f, g, nb"<< dom.Particles[i]->x(0)<<", "<<fx[i]<< ", "<<gx[i]<<", "<<nb[i]<<endl;
  }
  cout << endl<< "Derivatives"<<endl;
  for (int i = 0; i<dom.Particles.Size();i++) {
    double x = dom.Particles[i]->x(0);
    cout << "Analytical" << dom.Particles[i]->x(0)<<", "<<x*x<<endl;
    cout << dom.Particles[i]->x(0)<<", "<<dfx[i]<<endl;
  }
  
  return 0;
}
MECHSYS_CATCH
