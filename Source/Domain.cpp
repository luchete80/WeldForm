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
#include <chrono>
//#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include <ctime> //Clock

#include <vector>

#define MIN_PS_FOR_NBSEARCH		1.e-6//TODO: MOVE TO CLASS MEMBER

#include <set>

//https://stackoverflow.com/questions/19240540/dynamically-allocating-array-explain/19240932#19240932
template <typename T>
void Initiate (T ***mat,int row){
	//initialize 10 x 20 array:
	*mat = new T*[row];
	for (int i = 0; i < row; i++)
    (*mat)[i] = new T	; //DO NOT FORGET PARENTHESES, PRECEDENCE OF [] is HIGHER THAN DEREF *
	
	// for (int i=0;i<row;i++)
		// cout<< "row size"<<sizeof(*mat[i])/sizeof(float)<<endl;
}

using namespace std;

namespace SPH {
void General(Domain & dom)
{
}

void OutPut(Particle * Particles, double & Prop1, double & Prop2,  double & Prop3)
{
	Prop1 = 0.0;
	Prop2 = 0.0;
	Prop3 = 0.0;
}

void InFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, Boundary & bdry)
{
	Vel = bdry.inv;
	Den = bdry.inDensity;
}

void OutFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, Boundary & bdry)
{
	Vel = bdry.outv;
	Den = bdry.outDensity;
}

void AllFlowCon(Vec3_t & position, Vec3_t & Vel, double & Den, Boundary & bdry)
{
	Vel = bdry.allv;
	Den = bdry.allDensity;
}

// Constructor
inline Domain::Domain ()
{
    OutputName[0] = "Property1";
    OutputName[1] = "Property2";
    OutputName[2] = "Property3";
    Time    = 0.0;

    Dimension = 2;
    DomSize	= 0.0,0.0,0.0;

    Gravity	= 0.0,0.0,0.0;

    Cellfac = 2.0;

    KernelType	= 0;
    VisEq	= 0;
    Scheme	= 0;
		GradientType = 0;

    XSPH	= 0.0;
    InitialDist = 0.0;

    AvgVelocity = 0.0;
    hmax	= 0.0;

    omp_init_lock (&dom_lock);
    Nproc	= 1;

    deltat	= 0.0;
    deltatint	= 0.0;
    deltatmin	= 0.0;
    sqrt_h_a = 0.0025;

    TRPR = 0.0;
    BLPF = 0.0;

    InCon = & InFlowCon;
    OutCon = & OutFlowCon;
    AllCon = & AllFlowCon;
    GeneralBefore = & General;
    GeneralAfter = & General;
    UserOutput = & OutPut;

    DomMax = -100000000000.0;
    DomMin = 100000000000.0;
    I = OrthoSys::I;
	
	Vol=0.;
		auto_ts = true;
		
	
	gradKernelCorr = false;
	contact = false;
	contact_force_factor =1.;
	friction = 0.0;
	update_contact_surface = true;
	ts_nb_inc = 5;
  fric_type = Fr_Dyn;
  m_contact_forces_time = 0.; //TODO: MOVE TO ANOTHER CLASS
  m_forces_artifvisc_time = 0.;
  m_forces_momentum_time = 0.;
  m_forces_tensors_time = 0.;
  m_forces_update_time = 0.;
  m_scalar_prop = 0.;
  
  kin_energy_sum = int_energy_sum = 0.;
	
	thermal_solver = false;
  contact_mesh_auto_update = true;
  meshcount = 0;
}

inline Domain::~Domain ()
{
	size_t Max = Particles.Size();
	for (size_t i=1; i<=Max; i++)  Particles.DelItem(Max-i);
}

inline void Domain::Periodic_X_Correction(Vec3_t & x, double const & h, Particle * P1, Particle * P2)
{
	if (DomSize(0)>0.0) {if (x(0)>2*Cellfac*h || x(0)<-2*Cellfac*h) {(P1->CC[0]>P2->CC[0]) ? x(0) -= DomSize(0) : x(0) += DomSize(0);}}
	if (DomSize(1)>0.0) {if (x(1)>2*Cellfac*h || x(1)<-2*Cellfac*h) {(P1->CC[1]>P2->CC[1]) ? x(1) -= DomSize(1) : x(1) += DomSize(1);}}
	if (DomSize(2)>0.0) {if (x(2)>2*Cellfac*h || x(2)<-2*Cellfac*h) {(P1->CC[2]>P2->CC[2]) ? x(2) -= DomSize(2) : x(2) += DomSize(2);}}
}

inline void Domain::Kernel_Set(Kernels_Type const & KT)
{
	KernelType = KT;
	if (KernelType==2) Cellfac = 3.0; else Cellfac = 2.0;
}

inline void Domain::Viscosity_Eq_Set(Viscosity_Eq_Type const & VQ)
{
	VisEq = VQ;
}

inline void Domain::Gradient_Approach_Set(Gradient_Type const & GT)
{
	GradientType = GT;
}

inline void Domain::AdaptiveTimeStep()
{
	if (deltatint>deltatmin)
	{
		if (deltat<deltatmin)
			deltat		= 2.0*deltat*deltatmin/(deltat+deltatmin);
		else
			deltat		= deltatmin;
	}
	else
	{
		if (deltatint!=deltat)
			deltat		= 2.0*deltat*deltatint/(deltat+deltatint);
		else
			deltat		= deltatint;
	}
	
	// if (contact){
		// if (min_force_ts < deltat)
		// //cout << "Step size changed minimum Contact Forcess time: " << 	min_force_ts<<endl;
		// deltat = min_force_ts;
	// }

	// if (deltat<(deltatint/1.0e5))
		// //cout << "WARNING: Too small time step, please choose a smaller time step initially to make the simulation more stable"<<endl;
		// throw new Fatal("Too small time step, please choose a smaller time step initially to make the simulation more stable");
}

inline void Domain::CheckMinTSVel() {
  //Min time step check based on velocity
  double test	= 0.0;

  deltatmin	= deltatint;
  #pragma omp parallel for schedule (static) private(test) num_threads(Nproc)
  for (int i=0; i<Particles.Size(); i++) {
    if (Particles[i]->IsFree) {
      test = 0.4 * Particles[i]->h/(Particles[i]->Cs + norm(Particles[i]->v));
      if (deltatmin > test ) {
        omp_set_lock(&dom_lock);
          deltatmin = test;
        omp_unset_lock(&dom_lock);
      }
    }
  }
}

inline void Domain::AddSingleParticle(int tag, Vec3_t const & x, double Mass, double Density, double h, bool Fixed)
{
   	Particles.Push(new Particle(tag,x,Vec3_t(0,0,0),Mass,Density,h,Fixed));
}

inline void Domain::AddBoxLength(int tag, Vec3_t const & V, double Lx, double Ly, double Lz, 
									double r, double Density, double h, int type, int rotation, bool random, bool Fixed) {
    if ( !(type == 0 || type == 1) ) {
	   	std::cout << "Packing Type is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => Hexagonal Close Packing" << std::endl;
		std::cout << "1 => Cubic Packing" << std::endl;
	    abort();
    }

    if (!(rotation==0 || rotation==90)) {
	   	std::cout << "Packing Rotation Angle is out of range. Please correct it and run again" << std::endl;
		std::cout << "0 => " << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << " 0 0 0 0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << " 0 0 0 0" << std::endl;
		std::cout << std::endl;
		std::cout << "90 => Cubic Close Packing" << std::endl;
		std::cout << "  0   0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << "0 0 0 0" << std::endl;
		std::cout << "0   0  " << std::endl;
		abort();
    }

//	Util::Stopwatch stopwatch;
    std::cout << "\n--------------Generating particles by AddBoxLength with defined length of particles-----------" << std::endl;

    size_t PrePS = Particles.Size();

    double x,y,xp,yp;
    size_t i,j;

    double qin = 0.03;
    srand(100);
		
		//For new SOA accessing
		std::vector <Vec3_t> x_sta;

    if (Dimension==3) {
    	if (type==0) {
    		//Hexagonal close packing
    		double z,zp;
			size_t k=0;
			zp = V(2);

			while ( zp <= (V(2)+Lz-r) ) {
				
				j = 0;
				yp = V(1);
				while (yp <= (V(1)+Ly-r)) {
					i = 0;
					xp = V(0);
					while (xp <= (V(0)+Lx-r))
					{
						if ((k%2!=0) && (j%2!=0)) x = V(0) + (2*i+(j%2)+(k%2)-1)*r; else x = V(0) + (2*i+(j%2)+(k%2)+1)*r;
						y = V(1) + (sqrt(3.0)*(j+(1.0/3.0)*(k%2))+1)*r;
						z = V(2) + ((2*sqrt(6.0)/3)*k+1)*r;
						if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
						else    		{Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),0.0,Density,h,Fixed));
													x_sta.push_back(Vec3_t(x,y,z));
						}
						i++;
						if ((k%2!=0) && (j%2!=0)) xp = V(0) + (2*i+(j%2)+(k%2)-1)*r; else xp = V(0) + (2*i+(j%2)+(k%2)+1)*r;
					}
					j++;
					yp = V(1) + (sqrt(3.0)*(j+(1.0/3.0)*(k%2))+1)*r;
				}
				k++;
				zp = V(2) + ((2*sqrt(6.0)/3)*k+1)*r;
				//cout << "Z: "<<z<<endl;
			}
    	}
    	else {
    		//Cubic packing
    		double z,zp;
			size_t k=0;
			zp = V(2);

			while (zp <= (V(2)+Lz-r)) {
				j = 0;
				yp = V(1);
				while (yp <= (V(1)+Ly-r))
				{
					//cout << "Y: "<<yp<<endl;
					i = 0;
					xp = V(0);
					while (xp <= (V(0)+Lx-r))
					{
						x = V(0) + (2.0*i+1)*r;
						y = V(1) + (2.0*j+1)*r;
						z = V(2) + (2.0*k+1)*r;
						if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
						else    		Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),0.0,Density,h,Fixed));
						x_sta.push_back(Vec3_t(x,y,z));
						i++;
						xp = V(0) + (2*i+1)*r; //COMMENTED BY LUCIANO
						//cout << "X: "<<xp<<endl;
					}
					j++;
					yp = V(1) + (2.0*j+1)*r;//COMMENTED BY LUCIANO
				}
				k++;
				zp = V(2) + (2.0*k+1)*r;//COMMENTED BY LUCIANO
				cout << "Z: "<<z<<endl;
			}
    	}

        //Calculate particles' mass in 3D
        Vec3_t temp, Max=V;
		for (size_t i=PrePS; i<Particles.Size(); i++) {
			if (Particles[i]->x(0) > Max(0)) Max(0) = Particles[i]->x(0);
			if (Particles[i]->x(1) > Max(1)) Max(1) = Particles[i]->x(1);
			if (Particles[i]->x(2) > Max(2)) Max(2) = Particles[i]->x(2);
		}
		Max +=r;
		temp = Max-V;
		cout << "BoxDimensions: "<<temp(0)<<", "<<temp(1)<<", "<<temp(2)<<", "<<endl;
		double Mass = temp(0)*temp(1)*temp(2)*Density/(Particles.Size()-PrePS);
		
		cout << "Particle mass: " << Mass <<endl;
		
		// New SOA members
		cout << "Allocating "<<endl;
		Initiate (&m_x,Particles.Size());
		Initiate (&m_h,Particles.Size());
		Initiate (&m_kT,Particles.Size());
		Initiate (&m_cpT,Particles.Size());
		Initiate (&m_hcT,Particles.Size());
		Initiate (&m_qconvT,Particles.Size());
		Initiate (&m_T,Particles.Size());		
		Initiate (&m_Tinf,Particles.Size());
		Initiate (&m_dTdt,Particles.Size());
		Initiate (&m_rho,Particles.Size());
		Initiate (&m_mass,Particles.Size());
		cout << "Done."<<endl;
    sigma = new double [6*Particles.Size()];
    
		// m_x 		= new Vec3_t* [Particles.Size()];
		// m_h 		= new double* [Particles.Size()];
		// m_kT 		= new double* [Particles.Size()];
		// m_cpT 	= new double* [Particles.Size()];
		// m_hcT 	= new double* [Particles.Size()];
		// m_qconvT= new double* [Particles.Size()];
		// m_T 		= new double* [Particles.Size()];
		// m_Tinf 	= new double* [Particles.Size()];
		// m_dTdt 	= new double* [Particles.Size()];
		// m_rho 	= new double* [Particles.Size()];
		// m_mass 	= new double* [Particles.Size()];		
		//TODO-> CHANGE TO static members (particle will be deleted)
		for (int p=0;p<Particles.Size();p++){			
			m_x[p] 		= &x_sta[p];
			*m_rho[p] 	= Density; 
			*m_T[p] = *m_Tinf[p] = *m_hcT[p] = *m_kT[p] = *m_qconvT[p] = 0.;
			*m_mass[p] = Mass;
		}
		
		#pragma omp parallel for num_threads(Nproc)
		#ifdef __GNUC__
		for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
		#else
		for (int i=0; i<Particles.Size(); i++)//Like in Domain::Move
		#endif
		{
			Particles[i]->Mass = Mass;
		}
    } else if (Dimension==2) {
    	if (type==0)
    	{
    		//Hexagonal close packing
    		if (rotation==0)
    		{
				j = 0;
				yp = V(1);

				while (yp <= (V(1)+Ly-r))
				{
					i = 0;
					xp = V(0);
					while (xp <= (V(0)+Lx-r))
					{
						x = V(0) + (2*i+(j%2)+1)*r;
						y = V(1) + (sqrt(3.0)*j+1)*r;
						if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
							else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						i++;
						xp = V(0) + (2*i+(j%2)+1)*r;
					}
					j++;
					yp = V(1) + (sqrt(3.0)*j+1)*r;
				}
			}
    		else
    		{
				i = 0;
				xp = V(0);

				while (xp <= (V(0)+Lx-r))
				{
					j = 0;
					yp = V(1);
					while (yp <= (V(1)+Ly-r))
					{
						x = V(0) + (sqrt(3.0)*i+1)*r;
						y = V(1) + (2*j+(i%2)+1)*r;
						if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
							else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						j++;
						yp = V(1) + (2*j+(i%2)+1)*r;
					}
					i++;
					xp = V(0) + (sqrt(3.0)*i+1)*r;
				}
    		}
    	}
    	else
    	{
    		//Cubic packing
    		j = 0;
			yp = V(1);

			while (yp <= (V(1)+Ly-r))
			{
				i = 0;
				xp = V(0);
				while (xp <= (V(0)+Lx-r))
				{
					x = V(0) + (2*i+1)*r;
					y = V(1) + (2*j+1)*r;
					if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),0.0),Vec3_t(0,0,0),(sqrt(3.0)*r*r)*Density,Density,h,Fixed));
						else    Particles.Push(new Particle(tag,Vec3_t(x,y,0.0),Vec3_t(0,0,0),2.0*r*2.0*r*Density,Density,h,Fixed));
					i++;
					xp = V(0) + (2*i+1)*r;
				}
				j++;
				yp = V(1) + (2*j+1)*r;
			}

    	}
    }
		
		cout << "Particle Count: "<<Particles.Size()<< endl;

	R = r;
}

inline void Domain::Add3DCubicBoxParticles(int tag, Vec3_t const & V, double Lx, double Ly, double Lz, 
									double r, double Density, double h) {
//	Util::Stopwatch stopwatch;
    std::cout << "\n--------------Generating particles by AddBoxLength with defined length of particles-----------" << std::endl;

    double x,y,xp,yp;
    size_t i,j;
	size_t PrePS = Particles.Size();
    
	if (Dimension==3) {
		//Cubic packing
		double z,zp;
		size_t k=0;
		zp = V(2);

		while (zp <= (V(2)+Lz-r)) {
			
			j = 0;
			yp = V(1);
			while (yp <= (V(1)+Ly-r))
			{
				//cout << "Y: "<<yp<<endl;
				i = 0;
				xp = V(0);
				while (xp <= (V(0)+Lx-r))
				{
					x = V(0) + (2.0*i+1)*r;
					y = V(1) + (2.0*j+1)*r;
					z = V(2) + (2.0*k+1)*r;
					Particles.Push(new Particle(tag,Vec3_t(x,y,z),Vec3_t(0,0,0),0.0,Density,h,false));
					i++;
					xp = V(0) + (2*i+1)*r; //COMMENTED BY LUCIANO
					//cout << "X: "<<xp<<endl;
				}
				j++;
				yp = V(1) + (2.0*j+1)*r;//COMMENTED BY LUCIANO
			}
			k++;
			cout << "Z: "<<z<<endl;
			zp = V(2) + (2.0*k+1)*r;//COMMENTED BY LUCIANO
		}
    	//Vol+=(Lx*Ly*Lz);
        Vec3_t temp, Max=V;
		for (size_t i=PrePS; i<Particles.Size(); i++) {
			if (Particles[i]->x(0) > Max(0)) Max(0) = Particles[i]->x(0);
			if (Particles[i]->x(1) > Max(1)) Max(1) = Particles[i]->x(1);
			if (Particles[i]->x(2) > Max(2)) Max(2) = Particles[i]->x(2);
		}
		Max +=r;
		temp = Max-V;
		cout << "BoxDimensions: "<<temp(0)<<", "<<temp(1)<<", "<<temp(2)<<", "<<endl;
		Vol+=temp(0)*temp(1)*temp(2);
		//double Mass = temp(0)*temp(1)*temp(2)
    }//Dimension
}

// Calculate Mass for 3D particles
inline void Domain::Calculate3DMass(double Density){
	double Mass = Vol*Density/Particles.Size();
	cout << "Particle Mass: "<<Mass<<endl;
	#pragma omp parallel for num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
	#else
	for (int i=0; i<Particles.Size(); i++)//Like in Domain::Move
	#endif
	{
		Particles[i]->Mass = Mass;
	}

}
	
//////Return half (on the quadrant) particle count from a single position in an axis
int calcHalfPartCount(const double &r, const double &R, const int xinc){
	int ypartcount = -1;
	if ( xinc > 0 ){
		ypartcount = 1;
		double yp = r;
		double xp = r + (double)(xinc - 1 ) *2.*r; 
		double rad = sqrt(yp*yp + xp*xp);
		while( rad <= R -r ){
			yp += 2.*r;
			rad = sqrt(yp*yp + xp*xp);
			ypartcount++;
		}
		ypartcount-=1;
	}
	return ypartcount;
}

int Domain::AssignZone(Vec3_t &start, Vec3_t &end, int &id){
  int partcount = 0;
  for (size_t a=0; a<Particles.Size(); a++){
    bool included=true;
    for (int i=0;i<3;i++){
      if (Particles[a]->x(i) < start[i] || Particles[a]->x(i) > end[i])
        included = false;
    }
    if (included){
      Particles[a]->ID=id; 
      partcount++;
    }
  }
  return partcount;
}

inline void Domain::AddCylinderLength(int tag, Vec3_t const & V, double Rxy, double Lz, 
									double r, double Density, double h, bool Fixed, bool ghost) { //ghost refers to symmetry at bottom z coordinate

//	Util::Stopwatch stopwatch;
    std::cout << "\n--------------Generating particles by CylinderBoxLength with defined length of particles-----------" << std::endl;

    size_t PrePS = Particles.Size();

    double xp,yp;
    size_t i,j;

    double qin = 0.03;
    srand(100);
	
	double Lx, Ly;
	
	//Particles are tried to be aligned 
	int numpartxy=1;	//MAX ON EACH EDGE
	//PARTCILES ARE NOT ALIGNED WITH AXIS; BUT SYMMETRIC ON EACH QUADRANT
	//MIN CONFIG IS 4 PARTICLES; ALWAYS NUM PARTICLES IS PAIR
	numpartxy = calcHalfPartCount(r, Rxy, 1);
	
	//// GHOST THING

	int ghost_rows = 3; 

	int xy_ghost_part_count[ghost_rows];
	//cout << "X/Y Particles: " << numpartxy<<endl;
	//yp=pos;
	int numypart,numxpart;
	int xinc,yinc,yinc_sign;
	
	int id_part=0;
	
  if (Dimension==3) {
		int part_per_row = 0;
    	//Cubic packing
		double zp;
		size_t k=0;
		zp = V(2)/*+r*/;
		//Calculate row count for non ghost particles
		while (zp <= (V(2)+Lz -r)){
			k++; 
      zp += 2.*r;      
		}
		cout << "Particle Row count: "<< k << endl;
		int last_nonghostrow = k;
		
		k = 0;zp = V(2)/*+r*/;

		while (zp <= (V(2)+Lz -r)) {
			j = 0;
			yp = V(1) - r - (2.*r*(numpartxy - 1) ); //First increment is radius, following ones are 2r
			//cout << "y Extreme: "<<yp<<endl;
			
			numypart = 2*numpartxy;	//And then diminish by 2 on each y increment
			yinc = numpartxy;	//particle row from the axis
			yinc_sign=-1;
			//cout << "y max particles: "<<numypart<<endl;
			for (j=0;j<numypart;j++){
				//cout << "y inc: "<<yinc<<endl;
				numxpart = calcHalfPartCount(r, Rxy, yinc);
				//cout << "xpart: "<< numxpart<<endl;
				xp = V(0) - r - (2.*r*(numxpart - 1) ); //First increment is radius, following ones are 2r
				for (i=0; i<2*numxpart;i++) {
					//if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
					//	else    
					Particles.Push(new Particle(tag,Vec3_t(xp,yp,zp),Vec3_t(0,0,0),0.0,Density,h,Fixed));
					if (zp == V(2))
						part_per_row++;
					id_part++;
					xp += 2.*r;
				}
				yp += 2.*r;
				yinc+=yinc_sign;
				if (yinc<1) {//Reach the axis, now positive increments
					yinc = 1;
					yinc_sign=1;
				}
			}
			k++;
      zp += 2.*r;
		}
		cout << "Particles per row: "<<part_per_row<<endl;
		
    zp = V(2) - 2.0*r;
		//Insert ghost pairs relation
		if (ghost){
      //// Z PLANE, BOTTOM COORDINATE /////
      cout << "inserting z ghost particles at z bottom..."<<endl;
      //z Symm particles
      int sym_part;
      int part = 0;
      id_part = Particles.Size(); //REDUNDANT

      for (int zinc = 0; zinc < ghost_rows ;zinc++){
        for (int xy = 0; xy < part_per_row;xy++){
          xp = Particles[part]->x(0);
          yp = Particles[part]->x(1);
          //zp = Particles[part]->x(2);
          Particles.Push(new Particle(tag,Vec3_t(xp,yp,zp),Vec3_t(0,0,0),0.0,Density,h,Fixed));				
          Particles[id_part]->inner_mirr_part = part;
          //Particles[id_part]->ID = -50;
          //cout << "part , sym"<<part<<", "<<id_part<<endl;
          Particles[id_part]->ghost_plane_axis = 2;
          Particles[id_part]->not_write_surf_ID = true; //TO NOT BE WRITTEN BY OUTER SURFACE CALC
          
          //ONLY FOR TESTING SYMMETRY PLANES!
          //Particles[id_part]->ID = Particles[id_part]->ghost_plane_axis; //ONLY FOR TESTING IN PARAVIEW!   
          GhostPairs.Push(std::make_pair(part,id_part));
          
          id_part++;
          part++;
        }
        zp -= 2.0*r;
      }
		////// PARALLELIZE!
		}
		///////Calculate particles' mass in 3D
		// Vec3_t temp, Max=V;
		// for (size_t i=PrePS; i<Particles.Size(); i++) {
			// if (Particles[i]->x(0) > Max(0)) Max(0) = Particles[i]->x(0);
			// if (Particles[i]->x(1) > Max(1)) Max(1) = Particles[i]->x(1);
			// if (Particles[i]->x(2) > Max(2)) Max(2) = Particles[i]->x(2);
		// }
		// Max +=r;
		// temp = Max-V;
		// double Mass = temp(0)*temp(1)*temp(2)*Density/(Particles.Size()-PrePS);
		
		double Vol = M_PI * Rxy * Rxy * Lz;		
		//double Mass = Vol * Density / (Particles.Size()-PrePS);
		double Mass = Vol * Density /Particles.Size();
		
		cout << "Total Particle count: " << Particles.Size() <<endl;
		cout << "Particle mass: " << Mass <<endl;

		#pragma omp parallel for num_threads(Nproc)
		#ifdef __GNUC__
		for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
		#else
		for (int i=0; i<Particles.Size(); i++)//Like in Domain::Move
		#endif
		{
			Particles[i]->Mass = Mass;
		}

	}//Dim 3

	R = r;
}

void Domain::AddFixedMassScaling(const double &factor){
  mass_scaling_factor = factor;
 
  #pragma omp parallel for num_threads(Nproc)
  #ifdef __GNUC__
  for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
  #else
  for (int i=0; i<Particles.Size(); i++)//Like in Domain::Move
  #endif
  {
    Particles[i]->Mass *= factor;
    Particles[i]->Density *= factor;
  } 
}
//////////////////////////////////////
// HERE PARTICLE DISTRIBUTION IS RADIAL (DIFFERENT FROM PREVIOUS )
void Domain::AddQuarterCylinderLength(int tag, double Rxy, double Lz, 
																				double r, double Density, double h, bool Fixed, bool symlength = false) {
	//Util::Stopwatch stopwatch;
	std::cout << "\n--------------Generating particles by CylinderBoxLength with defined length of particles-----------" << std::endl;

	size_t PrePS = Particles.Size();
	double xp,yp;
	size_t i,j;
	double qin = 0.03;
	srand(100);
	
	double Lx, Ly;
	
	//Particles are tried to be aligned 
	int numpartxy = 1;	//MAX ON EACH EDGE
	//PARTCILES ARE NOT ALIGNED WITH AXIS; BUT SYMMETRIC ON EACH QUADRANT
	//MIN CONFIG IS 4 PARTICLES; ALWAYS NUM PARTICLES IS PAIR
	numpartxy = calcHalfPartCount(r, Rxy, 1);
	

	//yp=pos;
	int numypart,numxpart;
	int xinc,yinc;
	
	int id_part=0;
	int ghost_rows = 2;
	
	double z0;
	if (symlength) 	z0 = r;
	else						z0 = -Lz/2. - r; //CHECK: -Lz/2. - r or -Lz/2.?
	
	int part_per_row=0;
  std::vector <int> symm_x;
  std::vector <int> symm_y;
  
  if (Dimension==3) {
    	//Cubic packing
		double zp;
		size_t k=0;
		zp = z0;
		//Calculate row count for non ghost particles
		while (zp <= (z0+Lz -r)){
			k++; zp = z0 + (2.0*k+1)*r;			
		}
		//cout << "Particle Row count: "<< k << endl;
		int last_nonghostrow = k;
		k = 0;zp = z0;

		while (zp <= ( z0 + Lz - r)) {
			j = 0;
			//yp = - r - (2.*r*(numpartxy - 1) ); //First increment is radius, following ones are 2r
			yp = r; //First increment is radius, following ones are 2r
			//cout << "y Extreme: "<<yp<<endl;
			
			numypart = numpartxy;	//And then diminish by 2 on each y increment
			yinc = 1;	//particle row from the axis


			cout << "y max particles: "<<numypart<<endl;
			for (j=0;j<numypart;j++){
				//cout << "y inc: "<<yinc<<endl;
				numxpart = calcHalfPartCount(r, Rxy, yinc);
				//cout << "xpart: "<< numxpart<<endl;
				//xp = - r - (2.*r*(numxpart - 1) ); //First increment is radius, following ones are 2r
				//It is convenient to allocate now the ghost (symmetry) variables?
				xp = r; //First increment is radius, following ones are 2r
				for (i=0; i < numxpart;i++) {
					//if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
					//	else    
					Particles.Push(new Particle(tag,Vec3_t(xp,yp,zp),Vec3_t(0,0,0),0.0,Density,h,Fixed));
          if ( i < ghost_rows ){
            symm_x.push_back(id_part);
            //if (k==0) Particles[id_part]->ID = id_part; //ONLY FOR TESTING IN PARAVIEW!
          }
          if ( j < ghost_rows) {
            symm_y.push_back(id_part);
            //if (k==0) Particles[id_part]->ID = id_part; //ONLY FOR TESTING IN PARAVIEW!
					}
          if (zp == z0)
						part_per_row++;
					
					id_part++;
					xp += 2.*r;
				}
				yp += 2.*r;
				yinc +=1;

			}
			k++;
			zp = z0 + (2.0*k+1)*r;
		}
			cout << "Particles per row: "<<part_per_row<<endl;
    cout << " symmetric particles x, y :"<< symm_x.size()<<", "<<symm_y.size()<<endl;
		
		double Vol = M_PI * Rxy * Rxy * Lz;		
		//double Mass = Vol * Density / (Particles.Size()-PrePS);
		double Mass = Vol * Density /Particles.Size();
		
		cout << "Total Particle count: " << Particles.Size() <<endl;
		cout << "Particle mass: " << Mass <<endl;

		#pragma omp parallel for num_threads(Nproc)
		#ifdef __GNUC__
		for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
		#else
		for (int i=0; i<Particles.Size(); i++)//Like in Domain::Move
		#endif
		{
			Particles[i]->Mass = Mass;
		}

	}//Dim 3

	R = r;									
}


//////////////////////////////////////
// HERE PARTICLE DISTRIBUTION IS RADIAL (DIFFERENT FROM PREVIOUS )
void Domain::AddDoubleSymCylinderLength(int tag, double Rxy, double Lz, 
																				double r, double Density, double h, bool Fixed, bool symlength = false) {
	//Util::Stopwatch stopwatch;
	std::cout << "\n--------------Generating particles by CylinderBoxLength with defined length of particles-----------" << std::endl;

	size_t PrePS = Particles.Size();
	double xp,yp;
	size_t i,j;
	double qin = 0.03;
	srand(100);
	
	double Lx, Ly;
	
	//Particles are tried to be aligned 
	int numpartxy = 1;	//MAX ON EACH EDGE
	//PARTCILES ARE NOT ALIGNED WITH AXIS; BUT SYMMETRIC ON EACH QUADRANT
	//MIN CONFIG IS 4 PARTICLES; ALWAYS NUM PARTICLES IS PAIR
	numpartxy = calcHalfPartCount(r, Rxy, 1);
	

	//yp=pos;
	int numypart,numxpart;
	int xinc,yinc;
	
	int id_part=0;
	int ghost_rows = 2;
	
	double z0;
	if (symlength) 	z0 = r;
	else						z0 = -Lz/2. - r; //CHECK: -Lz/2. - r or -Lz/2.?
	
	int part_per_row=0;
  std::vector <int> symm_x;
  std::vector <int> symm_y;
  
  if (Dimension==3) {
    	//Cubic packing
		double zp;
		size_t k=0;
		zp = z0;
		//Calculate row count for non ghost particles
		while (zp <= (z0+Lz -r)){
			k++; zp += 2.0 * r;			
		}
		//cout << "Particle Row count: "<< k << endl;
		int last_nonghostrow = k;
		k = 0;zp = z0;

		while (zp <= ( z0 + Lz - r)) {
			j = 0;
			//yp = - r - (2.*r*(numpartxy - 1) ); //First increment is radius, following ones are 2r
			yp = r; //First increment is radius, following ones are 2r
			//cout << "y Extreme: "<<yp<<endl;
			
			numypart = numpartxy;	//And then diminish by 2 on each y increment
			yinc = 1;	//particle row from the axis


			cout << "y max particles: "<<numypart<<endl;
			for (j=0;j<numypart;j++){
				//cout << "y inc: "<<yinc<<endl;
				numxpart = calcHalfPartCount(r, Rxy, yinc);
				//cout << "xpart: "<< numxpart<<endl;
				//xp = - r - (2.*r*(numxpart - 1) ); //First increment is radius, following ones are 2r
				//It is convenient to allocate now the ghost (symmetry) variables?
				xp = r; //First increment is radius, following ones are 2r
				for (i=0; i < numxpart;i++) {
					//if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
					//	else    
					Particles.Push(new Particle(tag,Vec3_t(xp,yp,zp),Vec3_t(0,0,0),0.0,Density,h,Fixed));
          if ( i < ghost_rows ){
            symm_x.push_back(id_part);
            //if (k==0) Particles[id_part]->ID = id_part; //ONLY FOR TESTING IN PARAVIEW!
          }
          if ( j < ghost_rows) {
            symm_y.push_back(id_part);
            //if (k==0) Particles[id_part]->ID = id_part; //ONLY FOR TESTING IN PARAVIEW!
					}
          if (zp == z0)
						part_per_row++;
					
					id_part++;
					xp += 2.*r;
				}
				yp += 2.*r;
				yinc +=1;

			}
			k++;
			zp += 2.0 * r;
		}
			cout << "Particles per row: "<<part_per_row<<endl;
    cout << " symmetric particles x, y :"<< symm_x.size()<<", "<<symm_y.size()<<endl;
		//TODO: CONVERT THIS IN A QUARTER CYL SECTOR FUNCTION 
		//Is it convenient to allocate these particles at the end? 		
		//Allocate Symmetry particles, begining from x, y and z
		cout << "Creating ghost particles"<<endl;
    
    
    ///// X AND Y PLANE //////////
		zp = z0; k= 0;
		//cout << "zmax"<<( z0 + Lz - r)<<endl;
		while (zp <= ( z0 + Lz - r)) {

			numypart = numpartxy;	//And then diminish by 2 on each y increment
			yinc = 1;	//particle row from the axis

      int sym_y_count = 0;
      int sym_x_count;
			//cout << "y particles: "<<numypart<<endl;
			for (j=0; j < ghost_rows ; j++){
				xp = r;
				yp = - r - 2*r*(yinc -1); //First increment is radius, following ones are 2r			
				numxpart = calcHalfPartCount(r, Rxy, yinc);
				//cout << "x particles: "<<numypart<<endl;
        sym_x_count = j;
				for (i=0; i < numxpart;i++) {
					//if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
					//	else    
					Particles.Push(new Particle(tag,Vec3_t(xp,yp,zp),Vec3_t(0,0,0),0.0,Density,h,Fixed)); //First insert on y plane 
					Particles.Push(new Particle(tag,Vec3_t(yp,xp,zp),Vec3_t(0,0,0),0.0,Density,h,Fixed)); //Transpose for x plane

          //if (k==0) Particles[id_part  ]->ID = symm_y[sym_y_count]; //ONLY FOR TESTING IN PARAVIEW!  
          //SYMMETRY ON Y PLANE IS ALTERNATED ON INDICES
          //if (k==0) Particles[id_part+1]->ID = symm_x[sym_x_count]; //ONLY FOR TESTING IN PARAVIEW!           
          
					Particles[id_part  ]->inner_mirr_part = symm_y[sym_y_count];
					Particles[id_part+1]->inner_mirr_part = symm_x[sym_x_count];
          
          GhostPairs.Push(std::make_pair(symm_y[sym_y_count],id_part  ));
          GhostPairs.Push(std::make_pair(symm_x[sym_x_count],id_part+1));
          
          Particles[id_part  ]->ghost_plane_axis = 1;
          Particles[id_part+1]->ghost_plane_axis = 0;
          Particles[id_part  ]->is_ghost = true;
          Particles[id_part+1]->is_ghost = true;
          
          Particles[id_part]->not_write_surf_ID = true; //TO NOT BE WRITTEN BY OUTER SURFACE CALC
          Particles[id_part+1]->not_write_surf_ID = true; //TO NOT BE WRITTEN BY OUTER SURFACE CALC
          //ONLY FOR TESTING SYMMETRY PLANES!
          //Particles[id_part  ]->ID = 1; //ONLY FOR TESTING IN PARAVIEW!  
          //Particles[id_part+1]->ID = 0; //ONLY FOR TESTING IN PARAVIEW!   
					
          id_part+=2;
					xp += 2.*r;
          sym_y_count++;
          sym_x_count+=ghost_rows;
				}//x rows
				yinc++;
			}//y rows
			k++;
			zp += 2.0 * r;	
			//cout << "zp " <<zp<<endl;
		}
		
    //// Z PLANE, BOTTOM COORDINATE /////
    cout << "inserting z particles"<<endl;
    zp = z0 - 2.0*r;
		//Insert ghost pairs relation
		if (symlength){
      //// Z PLANE, BOTTOM COORDINATE /////
      cout << "inserting z ghost particles at z bottom..."<<endl;
      //z Symm particles
      int sym_part;
      int part = 0;
      id_part = Particles.Size(); //REDUNDANT

      for (int zinc = 0; zinc < ghost_rows ;zinc++){
        for (int xy = 0; xy < part_per_row;xy++){
          xp = Particles[part]->x(0);
          yp = Particles[part]->x(1);
          //zp = Particles[part]->x(2);
          Particles.Push(new Particle(tag,Vec3_t(xp,yp,zp),Vec3_t(0,0,0),0.0,Density,h,Fixed));				
          Particles[id_part]->inner_mirr_part = part;
          //Particles[id_part]->ID = -50;
          cout << "part , sym"<<part<<", "<<id_part<<endl;
          Particles[id_part]->ghost_plane_axis = 2;
          Particles[id_part]->not_write_surf_ID = true; //TO NOT BE WRITTEN BY OUTER SURFACE CALC
          Particles[id_part  ]->is_ghost = true;
          //ONLY FOR TESTING SYMMETRY PLANES!
          //Particles[id_part]->ID = Particles[id_part]->ghost_plane_axis; //ONLY FOR TESTING IN PARAVIEW!   
          GhostPairs.Push(std::make_pair(part,id_part));
          
          id_part++;
          part++;
        }
        zp -= 2.0*r;
      }
		////// PARALLELIZE!
		}
    
		
		double Vol = M_PI * Rxy * Rxy * Lz;		
		//double Mass = Vol * Density / (Particles.Size()-PrePS);
		double Mass = Vol * Density /Particles.Size();
		
		cout << "Total Particle count: " << Particles.Size() <<endl;
		cout << "Particle mass: " << Mass <<endl;

		#pragma omp parallel for num_threads(Nproc)
		#ifdef __GNUC__
		for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
		#else
		for (int i=0; i<Particles.Size(); i++)//Like in Domain::Move
		#endif
		{
			Particles[i]->Mass = Mass;
		}

	}//Dim 3

	R = r;									
}

inline void Domain::MoveGhost(){

	for (int gp=0; gp<GhostPairs.Size(); gp++){
		int  i = GhostPairs[gp].first;
		int gi = GhostPairs[gp].second;
		
    //ASSUMING SYMMETRY
		//See normal direction, if it is vertical
    // tg axis is the same speed
		Particles[gi]-> v  = Particles[i]-> v;
		Particles[gi]-> va = Particles[i]-> va;
		Particles[gi]-> vb = Particles[i]-> vb;
    
    int axis = Particles[gi]-> ghost_plane_axis;
    
		Particles[gi]-> v[axis]  = - Particles[i]-> v[axis];
		Particles[gi]-> va[axis] = - Particles[i]-> va[axis];
		Particles[gi]-> vb[axis] = - Particles[i]-> vb[axis];

		Particles[gi]-> a = 0.; //TO NOT INFLUENCE TIME STEP
		
		// Particles[gi]-> v[axis] = 	-Particles[gi]-> v[axis]		
		// Particles[gi]-> va[axis] = 	-Particles[gi]-> va[axis];
		// Particles[gi]-> vb[axis] = - Particles[gi]-> vb[axis];
		
		//Position (xghost + xi )/2 = x wall 
	}
}


inline void Domain::AddTractionProbeLength(int tag, Vec3_t const & V, double Rxy, double Lz_side,
											double Lz_neckmin,double Lz_necktot,double Rxy_center,
											double r, double Density, double h, bool Fixed) {

//	Util::Stopwatch stopwatch;
    std::cout << "\n--------------Generating particles by CylinderBoxLength with defined length of particles-----------" << std::endl;

    size_t PrePS = Particles.Size();

    double xp,yp;
    size_t i,j;

    double qin = 0.03;
    srand(100);
	
	double Lx, Ly, Lz;
	
	Lz = 2. * Lz_side + Lz_necktot;
	
	//Particles are tried to be aligned 
	int numpartxy;	//MAX ON EACH EDGE
	//PARTCILES ARE NOT ALIGNED WITH AXIS; BUT SYMMETRIC ON EACH QUADRANT
	//MIN CONFIG IS 4 PARTICLES; ALWAYS NUM PARTICLES IS PAIR
	
	
	//cout << "X/Y Particles: " << numpartxy<<endl;
	//yp=pos;
	int numypart,numxpart;
	int xinc,yinc,yinc_sign;
	double z_radiusreduction = (Lz_necktot-Lz_neckmin)/2.;
	double tan = (Rxy - Rxy_center)/(z_radiusreduction);
	double R;
	double z1 = V(2) + Lz_side - r;
	double z2 = V(2) + Lz_side + z_radiusreduction - r;
	double z3 = V(2) + Lz_side + z_radiusreduction + Lz_neckmin - r;
	double z4 = V(2) + Lz_side + Lz_necktot - r;
	
	int part = 0;
    if (Dimension==3) {
    	//Cubic packing
		double zp;
		size_t k=0;
		zp = V(2);
		bool center;
		while (zp <= ( V(2) + Lz - r )) {
			center = false;
			if 		( zp <= z1 || zp >  z4)		R = Rxy;
			else if ( zp > 	z1 && zp <= z2 )	R = Rxy - (zp - z1) * tan; 
			else if ( zp >= z2 && zp < z3 )		{R = Rxy_center; center = true;}
			else if ( zp >= z3 && zp < z4 )		R = Rxy_center + (zp - z3) * tan;
			
			
			numpartxy = calcHalfPartCount(r, R, 1);
			yp = V(1) - r - (2.*r*(numpartxy - 1) ); //First increment is radius, following ones are 2r
			//cout << "y Extreme: "<<yp<<endl;			
			numypart = 2*numpartxy;	//And then diminish by 2 on each y increment
			yinc = numpartxy;	//particle row from the axis
			yinc_sign=-1;
			//cout << "y max particles: "<<numypart<<endl;
			for (j=0;j<numypart;j++){
				//cout << "y inc: "<<yinc<<endl;
				numxpart = calcHalfPartCount(r, R, yinc);
				//cout << "xpart: "<< numxpart<<endl;
				xp = V(0) - r - (2.*r*(numxpart - 1) ); //First increment is radius, following ones are 2r
				for (i=0; i<2*numxpart;i++) {
					//if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
					//	else   
				if (center){
					//cout << "Particle "<<part<<", "<< Vec3_t(xp,yp,zp) <<endl;
				}						
					Particles.Push(new Particle(tag,Vec3_t(xp,yp,zp),Vec3_t(0,0,0),0.0,Density,h,Fixed));
					part++;
					xp += 2.*r;
				}
				yp += 2.*r;
				yinc+=yinc_sign;
				if (yinc<1) {//Reach the axis, now positive increments
					yinc = 1;
					yinc_sign=1;
				}
			}
			k++;
			zp = V(2) + (2.0*k+1)*r;
		}

//
// double Rxy, double Lz_side,
											// double Lz_neckmin,double Lz_necktot,double Rxy_center,
											// double r, double Density, double h, bool Fixed) {
												
		//Calculate particles' mass in 3D
		// Vec3_t temp, Max=V;
		// for (size_t i=PrePS; i<Particles.Size(); i++) {
			// if (Particles[i]->x(0) > Max(0)) Max(0) = Particles[i]->x(0);
			// if (Particles[i]->x(1) > Max(1)) Max(1) = Particles[i]->x(1);
			// if (Particles[i]->x(2) > Max(2)) Max(2) = Particles[i]->x(2);
		// }
		// Max +=r;
		// temp = Max-V;
		// double Mass = temp(0)*temp(1)*temp(2)*Density/(Particles.Size()-PrePS);

		double L_cone = ( Lz_necktot - Lz_neckmin )/2.;
		double Vol = 2.0 * Lz_side * M_PI * Rxy * Rxy + 
								       Lz_neckmin * M_PI * Rxy_center * Rxy_center +
							   2.0 * L_cone * M_PI * (Rxy_center * Rxy_center + Rxy * Rxy + Rxy * Rxy_center ) / 3.0 ; //Cones
		double Mass = Vol * Density / (Particles.Size()-PrePS);
		
		cout << "Particle mass: " << Mass <<endl;
		
		#pragma omp parallel for num_threads(Nproc)
		#ifdef __GNUC__
		for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
		#else
		for (int i=0; i<Particles.Size(); i++)//Like in Domain::Move
		#endif
		{
			Particles[i]->Mass = Mass;
		}

	}//Dim 3

	R = r;
}

// Calculate Free Surface (for contact and convection)
void Domain::CalculateSurface(const int &id){
	id_free_surf = id;
	double mi,mj;
	Particle *P1,*P2;
	Vec3_t xij;
	//TODO: SAVE THIS AT THE BEGINING AND PARALLELIZE
	double totmass=0.;
	for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
		totmass += Particles[i]->Mass;
		
	totmass /= Particles.Size();;
	//cout << "Totmass" <<	totmass <<endl;
	
	int maxid;
	if (contact)
		maxid = first_fem_particle_idx[0];
	else 
		first_fem_particle_idx[0] = Particles.Size();
	

	for (size_t i=0; i < maxid; i++)	{//Like in Domain::Move
		Particles[i] -> normal = 0.;
		Particles[i] -> ID = Particles [i] -> ID_orig;
	}
	
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t k=0; k<Nproc;k++) 
	#else
	for (int k=0; k<Nproc;k++) 
	#endif	
	{
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

	}//Nproc
	
	//Calculate Particle Neighbours
	
	//TODO: Parallelize with lock

	int surf_part =0;
	for (size_t i=0; i < maxid; i++)	{//Like in Domain::Move
	
		Particles[i]->normal *= 1./totmass;
		
		if ( norm(Particles[i]->normal) >= 0.25 * Particles[i]->h && Particles[i]->Nb <= 46) {//3-114 Fraser {
			if (!Particles[i]->not_write_surf_ID)
      Particles[i]->ID = id;
			surf_part++;
		}
	}
	//cout << "Surface particles: " << surf_part<<endl;
}


inline void Domain::DelParticles (int const & Tags)
{
    Array<int> idxs; // indices to be deleted

	#pragma omp parallel for schedule(static) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
	#else
	for (int i=0; i<Particles.Size(); i++)//Like in Domain::Move
	#endif
    {
        if (Particles[i]->ID == Tags)
		{
			omp_set_lock(&dom_lock);
        	idxs.Push(i);
			omp_unset_lock(&dom_lock);
		}
    }
    if (idxs.Size()<1) throw new Fatal("Domain::DelParticles: Could not find any particles to delete");
    Particles.DelItems (idxs);

    std::cout << "\n" << "Particle(s) with Tag No. " << Tags << " has been deleted" << std::endl;
}


inline void Domain::StartAcceleration (Vec3_t const & a) {

	#pragma omp parallel for schedule(static) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
	#else
	for (int i=0; i<Particles.Size(); i++)//Like in Domain::Move
	#endif
	{
	    if (Particles[i]->IsFree){
			// Tensile Instability for all soil and solid particles
			if (Particles[i]->TI > 0.0)
        		{
				// XY plane must be used, It is very slow in 3D
				if (Dimension == 2)
				{
					double teta, Sigmaxx, Sigmayy, C, S;

					if ((Particles[i]->Sigma(0,0)-Particles[i]->Sigma(1,1))!=0.0)
						teta = 0.5*atan(2.0*Particles[i]->Sigma(0,1)/(Particles[i]->Sigma(0,0)-Particles[i]->Sigma(1,1)));
					else
						teta = M_PI/4.0;

					C = cos(teta);
					S = sin(teta);
					Sigmaxx = C*C*Particles[i]->Sigma(0,0) + 2.0*C*S*Particles[i]->Sigma(0,1) + S*S*Particles[i]->Sigma(1,1);
					Sigmayy = S*S*Particles[i]->Sigma(0,0) - 2.0*C*S*Particles[i]->Sigma(0,1) + C*C*Particles[i]->Sigma(1,1);
					if (Sigmaxx>0) Sigmaxx = -Particles[i]->TI * Sigmaxx/(Particles[i]->Density*Particles[i]->Density); else Sigmaxx = 0.0;
					if (Sigmayy>0) Sigmayy = -Particles[i]->TI * Sigmayy/(Particles[i]->Density*Particles[i]->Density); else Sigmayy = 0.0;
					Particles[i]->TIR(0,0) = C*C*Sigmaxx + S*S*Sigmayy;
					Particles[i]->TIR(1,1) = S*S*Sigmaxx + C*C*Sigmayy;
					Particles[i]->TIR(0,1) = Particles[i]->TIR(1,0) = S*C*(Sigmaxx-Sigmayy);
				}
				else
				{
					Mat3_t Vec,Val,VecT,temp;
					double pc_ti_inv_d2=Particles[i]->TI/(Particles[i]->Density*Particles[i]->Density);//Precompute some values
					Rotation(Particles[i]->Sigma,Vec,VecT,Val);
					//Before
					// if (Val(0,0)>0) Val(0,0) = -Particles[i]->TI * Val(0,0)/(Particles[i]->Density*Particles[i]->Density); else Val(0,0) = 0.0;
					// if (Val(1,1)>0) Val(1,1) = -Particles[i]->TI * Val(1,1)/(Particles[i]->Density*Particles[i]->Density); else Val(1,1) = 0.0;
					// if (Val(2,2)>0) Val(2,2) = -Particles[i]->TI * Val(2,2)/(Particles[i]->Density*Particles[i]->Density); else Val(2,2) = 0.0;
					if (Val(0,0)>0) Val(0,0) = -pc_ti_inv_d2 * Val(0,0); else Val(0,0) = 0.0;
					if (Val(1,1)>0) Val(1,1) = -pc_ti_inv_d2 * Val(1,1); else Val(1,1) = 0.0;
					if (Val(2,2)>0) Val(2,2) = -pc_ti_inv_d2 * Val(2,2); else Val(2,2) = 0.0;

					Mult(Vec,Val,temp);
					Mult(temp,VecT,Particles[i]->TIR);
				}
			}
	    	}
	    	else
	    	{
	       		// Reset the pressure and the induced velocity for solid boundaries
	    		Particles[i]->NSv = 0.0;
	    		Particles[i]->Pressure = 0.0;
	        	set_to_zero(Particles[i]->Sigma);
	        	set_to_zero(Particles[i]->ShearStress);
	    	}



		//Reset to zero for all particles
		Particles[i]->a		= a;
		Particles[i]->SatCheck	= false;
		Particles[i]->dDensity	= 0.0;
		Particles[i]->VXSPH	= 0.0;
		Particles[i]->ZWab	= 0.0;
		Particles[i]->SumDen	= 0.0;
		Particles[i]->SumKernel	= 0.0;
		if (Dimension == 2) Particles[i]->v(2) = 0.0;
		set_to_zero(Particles[i]->StrainRate);
		set_to_zero(Particles[i]->RotationRate);
		
	}
}

inline void Domain::PrimaryComputeAcceleration () {
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t k=0; k<Nproc;k++) 
	#else
	for (int k=0; k<Nproc;k++) 
	#endif
	
	{
		size_t P1,P2;
		Vec3_t xij;
		double h,K;
		// Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
		for (size_t a=0; a<FSMPairs[k].Size();a++) {
			P1	= FSMPairs[k][a].first;
			P2	= FSMPairs[k][a].second;
			xij	= Particles[P1]->x-Particles[P2]->x;
			h	= (Particles[P1]->h+Particles[P2]->h)/2.0;

			Periodic_X_Correction(xij, h, Particles[P1], Particles[P2]);

			K	= Kernel(Dimension, KernelType, norm(xij)/h, h);
			if ( !Particles[P1]->IsFree ) {
				omp_set_lock(&Particles[P1]->my_lock);
										Particles[P1]->SumKernel+= K;
					Particles[P1]->Pressure	+= Particles[P2]->Pressure * K + dot(Gravity,xij)*Particles[P2]->Density*K;
					Particles[P1]->Sigma 	 = Particles[P1]->Sigma + K * Particles[P2]->Sigma;
					if (Particles[P1]->NoSlip)		Particles[P1]->NSv 	+= Particles[P2]->v * K;
				omp_unset_lock(&Particles[P1]->my_lock);
			} else {
				omp_set_lock(&Particles[P2]->my_lock);
										Particles[P2]->SumKernel+= K;
					Particles[P2]->Pressure	+= Particles[P1]->Pressure * K + dot(Gravity,xij)*Particles[P1]->Density*K;
					Particles[P2]->Sigma	 = Particles[P2]->Sigma + K * Particles[P1]->Sigma;
					if (Particles[P2]->NoSlip)		Particles[P2]->NSv 	+= Particles[P1]->v * K;
				omp_unset_lock(&Particles[P2]->my_lock);
			}
		}
	}

	// Calculateing the finala value of the smoothed pressure, velocity and stress for fixed particles
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t i=0; i<FixedParticles.Size(); i++)
	#else
	for (int i=0; i<FixedParticles.Size(); i++)
	#endif
  { 
    //for (int m=0;m<meshcount;m++)
    if (Particles[FixedParticles[i]]-> ID != contact_surf_id[0])  //ADDED TO Prevent adding surface (rigid contact) particles
		if (Particles[FixedParticles[i]]->SumKernel!= 0.0) {
			size_t a = FixedParticles[i];
			Particles[a]->Pressure	= Particles[a]->Pressure/Particles[a]->SumKernel;
			Particles[a]->Sigma	= 1.0/Particles[a]->SumKernel*Particles[a]->Sigma;
			if (Particles[a]->NoSlip)	Particles[a]->NSv	= Particles[a]->NSv/Particles[a]->SumKernel;

			// Tensile Instability for fixed soil and solid particles
			if (Particles[a]->TI > 0.0)
			{
				// XY plane must be used, It is very slow in 3D
				if (Dimension == 2) {
					double teta, Sigmaxx, Sigmayy, C, S;
					if ((Particles[a]->Sigma(0,0)-Particles[a]->Sigma(1,1))!=0.0)
						teta = 0.5*atan(2.0*Particles[a]->Sigma(0,1)/(Particles[a]->Sigma(0,0)-Particles[a]->Sigma(1,1)));
					else
						teta = M_PI/4.0;

					C = cos(teta);
					S = sin(teta);
					Sigmaxx = C*C*Particles[a]->Sigma(0,0) + 2.0*C*S*Particles[a]->Sigma(0,1) + S*S*Particles[a]->Sigma(1,1);
					Sigmayy = S*S*Particles[a]->Sigma(0,0) - 2.0*C*S*Particles[a]->Sigma(0,1) + C*C*Particles[a]->Sigma(1,1);
					if (Sigmaxx>0) Sigmaxx = -Particles[a]->TI * Sigmaxx/(Particles[a]->Density*Particles[a]->Density); else Sigmaxx = 0.0;
					if (Sigmayy>0) Sigmayy = -Particles[a]->TI * Sigmayy/(Particles[a]->Density*Particles[a]->Density); else Sigmayy = 0.0;
					Particles[a]->TIR(0,0) = C*C*Sigmaxx + S*S*Sigmayy;
					Particles[a]->TIR(1,1) = S*S*Sigmaxx + C*C*Sigmayy;
					Particles[a]->TIR(0,1) = Particles[a]->TIR(1,0) = S*C*(Sigmaxx-Sigmayy);
				}
				else {
					Mat3_t Vec,Val,VecT,temp;
					Rotation(Particles[a]->Sigma,Vec,VecT,Val);
					double pc_ti_inv_d2=Particles[a]->TI/(Particles[a]->Density*Particles[a]->Density);//Precompute some values
					// if (Val(0,0)>0) Val(0,0) = -Particles[a]->TI * Val(0,0)/(Particles[a]->Density*Particles[a]->Density); else Val(0,0) = 0.0;
					// if (Val(1,1)>0) Val(1,1) = -Particles[a]->TI * Val(1,1)/(Particles[a]->Density*Particles[a]->Density); else Val(1,1) = 0.0;
					// if (Val(2,2)>0) Val(2,2) = -Particles[a]->TI * Val(2,2)/(Particles[a]->Density*Particles[a]->Density); else Val(2,2) = 0.0;
					if (Val(0,0)>0) Val(0,0) = -pc_ti_inv_d2 * Val(0,0); else Val(0,0) = 0.0;
					if (Val(1,1)>0) Val(1,1) = -pc_ti_inv_d2 * Val(1,1); else Val(1,1) = 0.0;
					if (Val(2,2)>0) Val(2,2) = -pc_ti_inv_d2 * Val(2,2); else Val(2,2) = 0.0;

					Mult(Vec,Val,temp);
					Mult(temp,VecT,Particles[a]->TIR);
				}
			}
		}

  } 
}

inline void Domain::LastComputeAcceleration ()
{
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (int k=0; k<Nproc;k++) {
		for (size_t i=0; i<SMPairs[k].Size();i++)
			CalcForce2233(Particles[SMPairs[k][i].first],Particles[SMPairs[k][i].second]);

		for (int i=0; i<FSMPairs[k].Size();i++)
			CalcForce2233(Particles[FSMPairs[k][i].first],Particles[FSMPairs[k][i].second]);
	}
	
  m_clock_begin = clock();
	// CONTACT FORCES
	if (contact) {
		CalcContactForces();
		
	}
  m_contact_forces_time += (double)(clock() - m_clock_begin) / CLOCKS_PER_SEC;
  
		//Min time step check based on the acceleration
		double test	= 0.0;
		double test1 = 1000.;
    double test2 = 1000.;
    
		deltatmin	= deltatint;
		#pragma omp parallel for schedule (static) private(test,test1,test2) num_threads(Nproc)
		for (int i=0; i<Particles.Size(); i++) {
			if (Particles[i]->IsFree) {
				//test = sqrt(Particles[i]->h/norm(Particles[i]->a));
				//test1 = 1000.;
        test1 = sqrt_h_a * sqrt(Particles[i]->h/norm(Particles[i]->a));
        //cout << "time step with a criteria"<< test1<<endl;
				//test = 0.1 * Particles[i]->h/(Particles[i]->Cs + norm(Particles[i]->v));
				//if (norm(Particles[i]->v) != 0.){
        test2 = 1000.;
        //test2 = 0.1 * Particles[i]->h/(Particles[i]->Cs + norm(Particles[i]->v));
          //cout << "time step with v criteria"<< test2<<endl;
        //} else
        test = std::min(test1,test2);
				//if (deltatmin > (sqrt_h_a*test)) {
					if (deltatmin > test ) {
					omp_set_lock(&dom_lock);
						//deltatmin = sqrt_h_a*test
						deltatmin = test;
					omp_unset_lock(&dom_lock);
				}
			}
		}
  //cout << "deltatmin "<<deltatmin<<endl;
}
//
// First term of Eqn 3-68 kirk Fraser
void Domain::CalcKinEnergyEqn(){
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
	#else
	for (int i=0; i<Particles.Size(); i++)//Like in Domain::Move
	#endif
	{
		Particles[i]->dkin_energy_dt = 0.;
	}

  size_t P1,P2;	
	double temp;
	#pragma omp parallel for schedule (static) private (P1,P2) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t k=0; k<Nproc;k++) 
	#else
	for (int k=0; k<Nproc;k++) 
	#endif
	{
		Vec3_t xij,vij;
		double mi,mj,di,dj,GK;
		double h,K;
		// Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
		for (size_t a=0; a<SMPairs[k].Size();a++) {
			P1	= SMPairs[k][a].first;
			P2	= SMPairs[k][a].second;
			xij	= Particles[P1]->x-Particles[P2]->x;
			vij	= Particles[P1]->v-Particles[P2]->v;
			mi = Particles[P1]->Mass; mj = Particles[P2]->Mass;
			di = Particles[P1]->Density; dj = Particles[P2]->Density;

			h	= (Particles[P1]->h + Particles[P2]->h)/2.0;
			GK	= GradKernel(Dimension, KernelType, norm(xij)/h, h);	
			
			temp = 0.5 * mj*(Particles[P1]->Pressure/(di*di)+Particles[P2]->Pressure/(dj*dj))*
							dot(vij,GK*xij);
			 
			Particles[P1]->dkin_energy_dt +=temp; 
			Particles[P2]->dkin_energy_dt -=temp;			 
		}
	}
  double inc = 0.;
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
	#else
	for (int i=0; i<Particles.Size(); i++)//Like in Domain::Move
	#endif
	{
    //omp_set_lock(&dom_lock);
    inc += Particles[i]->dkin_energy_dt;
    //omp_unset_lock(&dom_lock);		
	}	
  kin_energy_sum += inc * deltat;
}

void Domain::CalcIntEnergyEqn(){
  double inc = 0.;
  int max;
  if (!contact) max = Particles.Size();
  else          max = first_fem_particle_idx[0];
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t i=0; i<max; i++)	//Like in Domain::Move
	#else
	for (int i=0; i<max; i++)//Like in Domain::Move
	#endif
	{
    Particles[i]->CalcIntEnergyEqn();
    //omp_set_lock(&dom_lock);
    inc += Particles[i]->dint_energy_dt;
    //omp_unset_lock(&dom_lock);		
	}
  
	int_energy_sum += inc * deltat;
}

//New, for Bonet gradient correction
inline void Domain::CalcGradCorrMatrix () {
	double di=0.0,dj=0.0,mi=0.0,mj=0.0;
	
	std::vector < Mat3_t> temp(Particles.Size());
	Mat3_t m,mt[2];
	
	//cout << "Applying grad corr"<<endl;
	//#pragma omp parallel for schedule (static) num_threads(Nproc) //LUCIANO: THIS IS DONE SAME AS PrimaryComputeAcceleration
	for ( size_t k = 0; k < Nproc ; k++) {
		Particle *P1,*P2;
		Vec3_t xij;
		double h,GK;
		//cout << "SMPairs[k].Size()"<<SMPairs[k].Size()<<endl;
		//TODO: DO THE LOCK PARALLEL THING
		for (size_t a=0; a<SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			//cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
			P1	= Particles[SMPairs[k][a].first];
			P2	= Particles[SMPairs[k][a].second];
			xij	= P1->x - P2->x;
			h	= (P1->h+P2->h)/2.0;
			GK	= GradKernel(Dimension, KernelType, norm(xij)/h, h);	
			
			di = P1->Density; mi = P1->Mass;
			dj = P2->Density; mj = P2->Mass;
		
			Dyad (Vec3_t(GK*xij),xij,m);
			mt[0] = mj/dj * m;
			mt[1] = mi/di * m;
			//cout << "mt"<<mt[0]<<endl;
			//omp_set_lock(&P1->my_lock);
			//SIGN IS NEGATIVE (IF POSITIVE, GRADIENT SIGN IS OPPOSITE)
			
			temp[SMPairs[k][a].first]  = temp[SMPairs[k][a].first]  - mt[0];  
			temp[SMPairs[k][a].second] = temp[SMPairs[k][a].second] - mt[1];
		}
	}//Nproc
	//cout << "Fixed Pairs"<<endl;
	for ( size_t k = 0; k < Nproc ; k++) {
		Particle *P1,*P2;
		Vec3_t xij;
		double h,GK;
		//TODO: DO THE LOCK PARALLEL THING
		//cout << "FSMPairs[k].Size()"<<FSMPairs[k].Size()<<endl;
		for (size_t a=0; a<FSMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			//cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
			P1	= Particles[FSMPairs[k][a].first];
			P2	= Particles[FSMPairs[k][a].second];
			xij	= P1->x - P2->x;
			h	= (P1->h+P2->h)/2.0;
			GK	= GradKernel(Dimension, KernelType, norm(xij)/h, h);	
			
			di = P1->Density; mi = P1->Mass;
			dj = P2->Density; mj = P2->Mass;
		
			Dyad (Vec3_t(GK*xij),xij,m);
			mt[0] = mj/dj * m;
			mt[1] = mi/di * m;
			//omp_set_lock(&P1->my_lock);
			//SIGN IS NEGATIVE (IF POSITIVE, GRADIENT SIGN IS OPPOSITE)
			
			temp[FSMPairs[k][a].first]  = temp[FSMPairs[k][a].first]  - mt[0];  
			temp[FSMPairs[k][a].second] = temp[FSMPairs[k][a].second] - mt[1];
		}
	}//Nproc
	//cout << "Inverting"<<endl;
	//#pragma omp parallel for schedule (static) num_threads(Nproc)	//LUCIANO//LIKE IN DOMAIN->MOVE
	//cout << "Inverting"<<endl;
	int max_id = Particles.Size();
	if (contact)
		max_id = first_fem_particle_idx[0];
		
	for (int i=0; i<max_id; i++){
		// cout << "part "<<i<<endl;
		//cout << "x: "<<Particles[i]->x<<endl;
		// cout << "nb: "<<Particles[i]->Nb<<endl;
		// if (!Particles[i]->IsFree) cout << "Fixed"<<endl;
		//cout << "temp "<<temp[i]<<endl;
		if (Dimension == 2)
			temp[i](2,2) = 1;
		/** Inverse.*/
		//inline void Inv (Mat3_t const & M, Mat3_t & Mi, double Tol=1.0e-10)}	
		if (Particles[i]->IsFree){
			Inv(temp[i],m);	

			Particles[i] ->gradCorrM = m;
			//cout << "Corr Matrix: " << m <<endl;
		} else {
			Particles[i] ->gradCorrM = I;
		}
	}	
}

//New, for Bonet gradient correction
inline void Domain::CalcGradCorrMixedMatrix () {
	double di=0.0,dj=0.0,mi=0.0,mj=0.0;
	
	std::vector < Mat3_t> temp(Particles.Size());
	Mat3_t m,mt[2];
	
	//cout << "Applying grad corr"<<endl;
	//#pragma omp parallel for schedule (static) num_threads(Nproc) //LUCIANO: THIS IS DONE SAME AS PrimaryComputeAcceleration
  double sumden[Particles.Size()];
  Vec3_t gamma [Particles.Size()];
  
	for ( size_t k = 0; k < Nproc ; k++) {
		Particle *P1,*P2;
		Vec3_t xij;
		double h,GK,K;
    int i,j;
		//cout << "SMPairs[k].Size()"<<SMPairs[k].Size()<<endl;
		//TODO: DO THE LOCK PARALLEL THING
		for (size_t a=0; a<SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			//cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
      
			P1	= Particles[SMPairs[k][a].first];
			P2	= Particles[SMPairs[k][a].second];
			xij	= P1->x - P2->x;
			h	= (P1->h+P2->h)/2.0;
			GK	= GradKernel(Dimension, KernelType, norm(xij)/h, h);				
      K	  =     Kernel(Dimension, KernelType, norm(xij)/h, h);
      
      sumden [i]+= mj/dj * K;
      sumden [j]-= mi/di * K;
      //Bonet Eqn 57  
      gamma [i] += GK*xij*mj/dj;
      gamma [j] -= GK*xij*mi/di;      
			
			di = P1->Density; mi = P1->Mass;
			dj = P2->Density; mj = P2->Mass;
		
			Dyad (Vec3_t(GK*xij),xij,m);
			mt[0] = mj/dj * m;
			mt[1] = mi/di * m;
			//cout << "mt"<<mt[0]<<endl;
			//omp_set_lock(&P1->my_lock);
			//SIGN IS NEGATIVE (IF POSITIVE, GRADIENT SIGN IS OPPOSITE)
			
			temp[SMPairs[k][a].first]  = temp[SMPairs[k][a].first]  - mt[0];  
			temp[SMPairs[k][a].second] = temp[SMPairs[k][a].second] - mt[1];
		}
	}//Nproc
	//cout << "Fixed Pairs"<<endl;
	for ( size_t k = 0; k < Nproc ; k++) {
		Particle *P1,*P2;
		Vec3_t xij;
		double h,GK;
		//TODO: DO THE LOCK PARALLEL THING
		//cout << "FSMPairs[k].Size()"<<FSMPairs[k].Size()<<endl;
		for (size_t a=0; a<FSMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			//cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
			P1	= Particles[FSMPairs[k][a].first];
			P2	= Particles[FSMPairs[k][a].second];
			xij	= P1->x - P2->x;
			h	= (P1->h+P2->h)/2.0;
			GK	= GradKernel(Dimension, KernelType, norm(xij)/h, h);	
			
			di = P1->Density; mi = P1->Mass;
			dj = P2->Density; mj = P2->Mass;
		
			Dyad (Vec3_t(GK*xij),xij,m);
			mt[0] = mj/dj * m;
			mt[1] = mi/di * m;
			//omp_set_lock(&P1->my_lock);
			//SIGN IS NEGATIVE (IF POSITIVE, GRADIENT SIGN IS OPPOSITE)
			
			temp[FSMPairs[k][a].first]  = temp[FSMPairs[k][a].first]  - mt[0];  
			temp[FSMPairs[k][a].second] = temp[FSMPairs[k][a].second] - mt[1];
		}
	}//Nproc
	//cout << "Inverting"<<endl;
	//#pragma omp parallel for schedule (static) num_threads(Nproc)	//LUCIANO//LIKE IN DOMAIN->MOVE
	//cout << "Inverting"<<endl;
	for (int i=0; i<Particles.Size(); i++){
		// cout << "part "<<i<<endl;
		//cout << "x: "<<Particles[i]->x<<endl;
		// cout << "nb: "<<Particles[i]->Nb<<endl;
		// if (!Particles[i]->IsFree) cout << "Fixed"<<endl;
		//cout << "temp "<<temp[i]<<endl;
		if (Dimension == 2)
			temp[i](2,2) = 1;
		/** Inverse.*/
		//inline void Inv (Mat3_t const & M, Mat3_t & Mi, double Tol=1.0e-10)}	
		if (Particles[i]->IsFree){
			Inv(temp[i],m);	

			Particles[i] ->gradCorrM = m;
			//cout << "Corr Matrix: " << m <<endl;
		} else {
			Particles[i] ->gradCorrM = I;
		}
	}	
}

inline void Domain::Move (double dt) {
	//cout << "BEGIN MOVE Time " << Time << "------------------------------"<<endl;
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (int i=0; i<Particles.Size(); i++)
		if (Particles[i]->IsFree) {
			if (Particles[i]->InOut>0) {
				Particles[i]->a = 0.0;
				if (Particles[i]->InOut == 1) {
					Particles[i]->dDensity = 0.0;
					Particles[i]->ZWab = 0.0;
				} else {
					if (BC.outDensity>0.0) {
						Particles[i]->dDensity = 0.0;
						Particles[i]->ZWab = 0.0;
					}
				}
			}
			//cout << "Particle: "<<i<<endl;
      //if (Particles[i]->is_ghost)
        Particles[i]->Move(dt,DomSize,TRPR,BLPF,Scheme,I);
      // if (i==624){
        // if (Particles[i]->eff_strain_rate>0)
          // cout << "particle 624, eff strain rate : "<<Particles[i]->eff_strain_rate<<", Et: "<<Particles[i]->Et<<"sigmaeq"<< Particles[i]->Sigma_eq<<", yield"<<Particles[i]->Sigmay<<endl;
      // }
		}
}

inline void Domain::WholeVelocity() {
    //Apply a constant velocity to all particles in the initial time step
    if (norm(BC.allv)>0.0 || BC.allDensity>0.0) {
    	Vec3_t vel = 0.0;
    	double den = 0.0;

	#pragma omp parallel for schedule (static) private(vel,den) num_threads(Nproc)
    	for (int i=0 ; i<Particles.Size() ; i++) {
		AllCon(Particles[i]->x,vel,den,BC);
    		if (Particles[i]->IsFree && norm(BC.allv)>0.0) {
			Particles[i]->v		= vel;
 		}
    		if (Particles[i]->IsFree && BC.allDensity>0.0) {
			Particles[i]->Density	= den;
			Particles[i]->Pressure	= EOS(Particles[i]->PresEq, Particles[i]->Cs, Particles[i]->P0,Particles[i]->Density, Particles[i]->RefDensity);
    		}
    	}
    }
}

inline void Domain::InitialChecks() {
	//initializing identity matrix
	if (Dimension == 2) I(2,2) = 0;

	if (Dimension<=1 || Dimension>3) {
		std::cout << "Please correct the dimension (2=>2D or 3=>3D) and run again" << std::endl;
		abort();
	}

	if (BC.InOutFlow>0 && BC.Periodic[0])
		throw new Fatal("Periodic BC in the X direction cannot be used with In/Out-Flow BC simultaneously");


	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for (int i=0; i<Particles.Size(); i++) //Initializing pressure of solid and fluid particles
			Particles[i]->Pressure = EOS(Particles[i]->PresEq, Particles[i]->Cs, Particles[i]->P0,Particles[i]->Density, Particles[i]->RefDensity);
}

inline void Domain::TimestepCheck ()
{
	// Check the time step
	double t1,t2;
	t1 = 0.25*hmax/(CsMax);
	if (MuMax>0.0) t2 = 0.125*hmax*hmax*rhomax/MuMax; else t2 =1000000.0;

	std::cout << "Max allowable time step using CFL = "<< std::min(t1,t2) << " S" << std::endl;
	std::cout << "User Time Step = "<< deltatint  << " S" << std::endl;

	// if (deltatint > std::min(t1,t2))
	// throw new Fatal("Please decrease the time step to the allowable range");
}

inline void Domain::Solve (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx) {
	std::cout << "\n--------------Solving---------------------------------------------------------------" << std::endl;

	size_t idx_out = 1;
	double tout = Time;

	//Initializing adaptive time step variables
	deltat = deltatint = deltatmin	= dt;
	
	auto start_whole = std::chrono::steady_clock::now();
		
	cout << "Initial checks"<<endl;
	InitialChecks();
	cout << "Creating list and cells"<<endl;
	CellInitiate();
	ListGenerate();
	PrintInput(TheFileKey);
	cout << "Time Step check"<<endl;
	TimestepCheck();
	WholeVelocity();
	
	//TODO: MOVE
	if (contact){
		for (int i=0; i<Particles.Size(); i++)
			Particles [i] -> ID_orig = Particles [i] -> ID;
	}
	
	cout << "Cell Size: "<<CellSize<<endl;
	
	if (contact) { //Calculate particle Stiffness
		//Cs	= sqrt(K/rho);
		for (int i=0; i<Particles.Size(); i++){
			double bulk = Particles[i]->Cs * Particles[i]->Cs *Particles[i]-> Density;  //RESTORE ORIGINAL BULK
			//TODO: If convection heat is updated every step, maybe dS could be calculated once 
			// in order to account for this too
			double dS = pow(Particles[i]->Mass/Particles[i]->Density,0.33333); //Fraser 3-119
			//Fraser Thesis, Eqn. 3-153
			Particles [i] -> cont_stiff = 9. * bulk * Particles [i]->G / (3. * bulk + Particles [i]->G) * dS; 
		}		
		cout << "dS, Cs, Contact Stiffness" << pow(Particles[0]->Mass/Particles[0]->Density,0.33333)<< ", " 
    << Particles[0]->Cs << ", " << Particles [0] -> cont_stiff <<endl;
		min_force_ts = deltat;
	}
	cout << "Fixed Particles Size: "<<FixedParticles.Size()<<endl;
	cout << "Initial Cell Number: "<<CellNo[0]<<", " <<CellNo[1]<<", "<< CellNo[2]<<", " <<endl;
	
	std::chrono::duration<double> total_time,neighbour_time;
	
	clock_t clock_beg;
	double clock_time_spent,start_acc_time_spent, pr_acc_time_spent,acc_time_spent, 
				contact_time_spent, trimesh_time_spent, bc_time_spent,
				mov_time_spent;

	double neigbour_time_spent_per_interval=0.;
	
	clock_time_spent = 
	pr_acc_time_spent=acc_time_spent= start_acc_time_spent = 
	contact_time_spent = trimesh_time_spent = bc_time_spent = 
	mov_time_spent = 0.;
  
  double contact_nb_time_spent = 0.;
  double contact_surf_time_spent = 0.;

	//Initial model output
	if (TheFileKey!=NULL) {
		String fn;
		fn.Printf    ("%s_Initial", TheFileKey);
		WriteXDMF    (fn.CStr());
		std::cout << "\nInitial Condition has been generated\n" << std::endl;
	}
	

	unsigned long steps=0;
	unsigned int first_step;
	
	int ts_i=0;

	bool isfirst = true;
	bool isyielding = false;

	//In case of contact this must be SURFACE particles
	//TODO, REMOVE so many nb search
  cout << "Calculating Nbs.."<<endl;
	if (contact){
		MainNeighbourSearch();
		SaveNeighbourData();				//Necesary to calulate surface! Using Particle->Nb (count), could be included in search
		CalculateSurface(1);				//After Nb search	
	}
  cout << "done."<<endl;
	//IF GRADCORR IS CALCULATED HERE; INVERSE IS NOT FOUND (ERROR)
	// TO BE CHECK
	// if (gradKernelCorr)
		// CalcGradCorrMatrix();	
	ClearNbData();
	
	//Print history
	std::ofstream of("History.csv", std::ios::out);
  of << "Displacement, pl_strain, eff_strain_rate, sigma_eq, sigmay, contforcesum"<<endl;
  
  bool check_nb_every_time = false;
  
  //TODO: Move to a function, AND HAS TO BE ZEROED EVERY NB SEARCH
  for(size_t i=0 ; i<Nproc ; i++) {
    std::vector<size_t> v(Particles.Size());
    ipair_SM.push_back(v); 
    jpair_SM.push_back(v); 
  }
  

	while (Time<=tf && idx_out<=maxidx) {
		clock_beg = clock();
		StartAcceleration(Gravity);
		start_acc_time_spent = (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		//if (BC.InOutFlow>0) InFlowBCFresh();
		auto start_task = std::chrono::system_clock::now();
		

		double max = 0;
		int imax;
		#pragma omp parallel for schedule (static) num_threads(Nproc)	//LUCIANO//LIKE IN DOMAIN->MOVE
		for (int i=0; i<Particles.Size(); i++){
			if (Particles[i]->pl_strain > max){
        omp_set_lock(&dom_lock);
				max= Particles[i]->pl_strain;
        omp_unset_lock(&dom_lock);
				imax=i;
			}
		}

    Vec3_t max_disp = Vec3_t(0.,0.,0.);
		for (int i=0; i<Particles.Size(); i++){
      for (int j=0;j<3;j++)
        if (Particles[i]->Displacement[j]>max_disp[j]){
          max_disp[j] = Particles[i]->Displacement [j];
          imax=i;
			}
		}
    
    // ATTENTION! COULD BE LARGE DISPLACEMENTS AND SMALL STRAINS 
    //EXAMPLE COMPRESSION WITH NO FRICTION, SO CONTACTS NBs SHOULD BE RECALCULATED
    if (norm(max_disp) > 0.1 * hmax){
      if (!check_nb_every_time)
        cout << "Checking Nb Every step now."<<endl;
      check_nb_every_time = true;
    }
    else 
      check_nb_every_time = false;
		
		if (max > MIN_PS_FOR_NBSEARCH && !isyielding){ //First time yielding, data has not been cleared from first search
			ClearNbData();
      
      // THIS IS IF MAINNBSEARCH INCLUDE SEARCHING CONTACT (NEW)
      // if (contact){
				// SaveNeighbourData();				//Necesary to calulate surface! Using Particle->Nb (count), could be included in search
				// CalculateSurface(1);				//After Nb search			        
      // }
      
			MainNeighbourSearch/*_Ext*/();
      
     // if (contact) SaveContNeighbourData();
			
			if (contact) {
				//TODO: CHANGE CONTACT STIFFNESS!
				SaveNeighbourData();				//Necesary to calulate surface! Using Particle->Nb (count), could be included in search
				CalculateSurface(1);				//After Nb search			
				ContactNbSearch();
				SaveContNeighbourData();	//Again Save Nb data
			}//contact
			isyielding  = true ;
		}
		if ( max > MIN_PS_FOR_NBSEARCH || isfirst || check_nb_every_time){	//TO MODIFY: CHANGE
			if ( ts_i == 0 ){
				clock_beg = clock();
				if (m_isNbDataCleared){

          // if (contact){
            // SaveNeighbourData();				//Necesary to calulate surface! Using Particle->Nb (count), could be included in search
            // CalculateSurface(1);				//After Nb search			        
          // }
					MainNeighbourSearch/*_Ext*/();
          //if (contact) SaveContNeighbourData();
					
          neigbour_time_spent_per_interval += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;					
          
          // TODO: SEPARATE CONTACT SEARCH STEP INTERVAL
					// OLD
          if (contact) {

						//cout << "performing contact search"<<endl
            clock_beg = clock();
          //if (update_contact_surface){
            
            SaveNeighbourData();				//Necesary to calulate surface! Using Particle->Nb (count), could be included in search
            CalculateSurface(1);				//After Nb search			
            contact_surf_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
            if (isfirst)
              CalcContactInitialGap(); //BEFORE! contactnb
            ContactNbSearch();
            SaveContNeighbourData();
            contact_nb_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
						//}
					}//contact				
				}// ts_i == 0				
				
			}
		
    } //( max > MIN_PS_FOR_NBSEARCH || isfirst ){	//TO MODIFY: CHANGE

		//NEW, gradient correction
			if (isfirst) {
				if (gradKernelCorr){
          cout << "Calculating gradient correction matrix"<<endl;
          CalcGradCorrMatrix();	}
				cout << "Done."<<endl;
				isfirst = false;
			}		

			
		auto end_task = std::chrono::system_clock::now();
		 neighbour_time = /*std::chrono::duration_cast<std::chrono::seconds>*/ (end_task- start_task);
		//std::cout << "neighbour_time (chrono, clock): " << clock_time_spent << ", " << neighbour_time.count()<<std::endl;
		clock_beg = clock();
		GeneralBefore(*this);
		bc_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		clock_beg = clock();
		PrimaryComputeAcceleration();
		pr_acc_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		clock_beg = clock();
		LastComputeAcceleration();
		acc_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		clock_beg = clock();

    double maxT = 0.;
    double minT =1000.;    
    for (size_t i=0; i<Particles.Size(); i++){
			//Particles[i]->T+= dt*Particles[i]->dTdt;
			Particles[i]->TempCalcLeapfrog(dt);
			if (Particles[i]->T > maxT)
				maxT=Particles[i]->T;
			if (Particles[i]->T < minT)
				minT=Particles[i]->T;      
    }
    

		if (thermal_solver){
			CalcConvHeat();
			CalcPlasticWorkHeat(deltat);
			CalcTempInc();
			CalcThermalExpStrainRate();	//Add Thermal expansion Strain Rate Term	
		}
		
		clock_beg = clock();
		GeneralAfter(*this);
		bc_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		steps++;
		//cout << "steps: "<<steps<<", time "<< Time<<", tout"<<tout<<endl;
		// output
		if (Time>=tout){
			if (TheFileKey!=NULL) {
				String fn;
				fn.Printf    ("%s_%04d", TheFileKey, idx_out);
				WriteXDMF    (fn.CStr());
				//fn.Printf    ("%s_%.5f", TheFileKey, Time);
				WriteCSV    (fn.CStr());

			}
			idx_out++;
			tout += dtOut;
			total_time = std::chrono::steady_clock::now() - start_whole;
			std::cout << "\nOutput No. " << idx_out << " at " << Time << " has been generated" << std::endl;
			std::cout << "Current Time Step = " <<deltat<<std::endl;
			
			clock_time_spent += neigbour_time_spent_per_interval;
			std::cout << "Total CPU time: "<<total_time.count() << endl <<
			", Nb: " << clock_time_spent << ", StAcc:  " 
			<< start_acc_time_spent << ", PrAcc: " << pr_acc_time_spent << endl<<
      "Total Forces : " << acc_time_spent<< endl<<", Contact Forces: "<< m_contact_forces_time<<
      "Artif Visc: "<<m_forces_artifvisc_time << ", Momentum forces: "<<m_forces_momentum_time << 
      "Forces Tensor: "<<m_forces_tensors_time<< endl<<
      " Forces Update: " << m_forces_update_time <<endl<<", Contact Nb : "<< contact_nb_time_spent << 
      "Contact Surf : "<< contact_surf_time_spent  << "Msh: " << trimesh_time_spent <<
			", BC: "<< bc_time_spent << 
			", mv: "<<mov_time_spent <<
      ", Contact Force Sum "<<contact_force_sum<<
      ", UserDefProp: "<<m_scalar_prop<<
			std::endl;
						
			cout << "Max plastic strain: " <<max<< "in particle" << imax << endl;
			
			std::cout << "Steps count in this interval: "<<steps-first_step<<"Total Step count"<<steps<<endl;
			cout << "Total Nb search time in this interval: " << neigbour_time_spent_per_interval;
			cout << "Average Nb search time in this interval: " << neigbour_time_spent_per_interval/(float)(steps-first_step)<<endl;

			cout << "Avg Neighbour Count"<<AvgNeighbourCount()<<endl;
			std::cout << "Max, Min, Avg temps: "<< maxT << ", " << minT << ", " << (maxT+minT)/2. <<std::endl;      
      cout << "Particle 0 pos and vel "<<endl;
      cout << Particles[0]->x<<endl;
      cout << Particles[0]->v<<endl;
      
      // cout << "ghost pair 0" << GhostPairs[0].first<<", "<<GhostPairs[0].second<<endl;
      // cout << Particles[GhostPairs[0].second]->x<<endl;
      // cout << Particles[GhostPairs[0].second]->v<<endl;
      
			first_step=steps;
			neigbour_time_spent_per_interval=0.;
			cout << "Max Displacements: "<<max_disp<<endl;
      
			if (contact)
				cout << "Max Contact Force: "<<max_contact_force<<endl;
			
			for (int p=0;p<Particles.Size();p++){
				if (Particles[p]->print_history)
          of << Particles[p]->Displacement << ", "<<Particles[p]->pl_strain<<", "<<Particles[p]->eff_strain_rate<<", "<< 
          Particles[p]->Sigma_eq<<", "  <<  Particles[p]->Sigmay << ", " <<
          contact_force_sum << endl;
			}
		}
		
		// if (isyielding)
			// cout << "Current Time Step: "<<deltat<<endl;
		
		// for (int i=0; i<Particles.Size(); i++){
			// if (Particles[i]->contforce>0.)
		if (auto_ts)
			AdaptiveTimeStep();
    //cout << "delta t"<<deltat<<endl;
    
		clock_beg = clock();

		Move(deltat); // INCLUDES GHOST PARTICLES
    MoveGhost();  //If Symmetry, 
    
    
		mov_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		clock_beg = clock();
		// Update velocity, plane coeff pplane and other things
    if (contact){
 		//cout << "checking contact"<<endl;
      if (contact_mesh_auto_update){
        for (int m=0; m<trimesh.size();m++)
          trimesh[m]->Update (deltat); //Update Node Pos, NOW includes PosCoeff and normals        
      }
      //cout << "Updating contact particles"<<endl;
      UpdateContactParticles(); //Updates normal and velocities
		}
    //cout << "Done"<<endl;

		trimesh_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		
		Time += deltat;
		//if (BC.InOutFlow>0) InFlowBCLeave(); else CheckParticleLeave ();
		
		
		if (max>MIN_PS_FOR_NBSEARCH){	//TODO: CHANGE TO FIND NEIGHBOURS
			if ( ts_i == (ts_nb_inc - 1) ){
				ClearNbData();
			}

			ts_i ++;
			if ( ts_i > (ts_nb_inc - 1) ) 
				ts_i = 0;
		
		}
		
	
	}
	

	of.close();
	
	std::cout << "\n--------------Solving is finished---------------------------------------------------" << std::endl;

}

void ContactNbUpdate(SPH::Domain *dom){
  dom->CalculateSurface(1);				//After Nb search			
  dom->ContactNbSearch();
  dom->SaveContNeighbourData();	//Again Save Nb data
}

// THIS IS LIKE THE FRASER ALGORITHM, LIKE STANDARD VERLET BUT X IS CALCULATED AFTER V
inline void Domain::SolveDiffUpdateModEuler (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx) {
	std::cout << "\n--------------Solving---------------------------------------------------------------" << std::endl;

	size_t idx_out = 1;
	double tout = Time;

	//Initializing adaptive time step variables
	deltat = deltatint = deltatmin	= dt;
	

	InitialChecks();
	CellInitiate();
	ListGenerate();
	PrintInput(TheFileKey);
	TimestepCheck();
	WholeVelocity();
	

	//Initial model output
	if (TheFileKey!=NULL) {
		String fn;
		fn.Printf    ("%s_Initial", TheFileKey);
		WriteXDMF    (fn.CStr());
		std::cout << "\nInitial Condition has been generated\n" << std::endl;
	}
	

	unsigned long steps=0;
	unsigned int first_step;
	
	int ts_i=0;

	bool isfirst = true;
	bool isyielding = false;
  
  cout << "Nb search"<<endl;
  ClearNbData();
  cout << "cleared"<<endl;
  MainNeighbourSearch();
  //SaveNeighbourData();				//Necesary to calulate surface! Using Particle->Nb (count), could be included in search
  cout << "Done"<<endl;
  //CalculateSurface(1);				//After Nb search	
	//ClearNbData();
	
	//Print history
	std::ofstream of("History.csv", std::ios::out);
  of << "Displacement, pl_strain, eff_strain_rate, sigma_eq, sigmay, contforcesum"<<endl;
  
  bool check_nb_every_time = false;

  cout << "Main Loop"<<endl;
  
	while (Time<=tf && idx_out<=maxidx) {
    
		StartAcceleration(Gravity);

		double max = 0;
		int imax;
		#pragma omp parallel for schedule (static) num_threads(Nproc)	//LUCIANO//LIKE IN DOMAIN->MOVE
		for (int i=0; i<Particles.Size(); i++){
			if (Particles[i]->pl_strain > max){
        omp_set_lock(&dom_lock);
				max= Particles[i]->pl_strain;
        omp_unset_lock(&dom_lock);
				imax=i;
			}
		}

    Vec3_t max_disp = Vec3_t(0.,0.,0.);
		for (int i=0; i<Particles.Size(); i++){
      for (int j=0;j<3;j++)
        if (Particles[i]->Displacement[j]>max_disp[j]){
          max_disp[j] = Particles[i]->Displacement [j];
          imax=i;
			}
		}
    
    // ATTENTION! COULD BE LARGE DISPLACEMENTS AND SMALL STRAINS 
    //EXAMPLE COMPRESSION WITH NO FRICTION, SO CONTACTS NBs SHOULD BE RECALCULATED
    // if (norm(max_disp) > 0.1 * hmax){
      // if (!check_nb_every_time)
        // cout << "Checking Nb Every step now."<<endl;
      // check_nb_every_time = true;
    // }
    // else 
      // check_nb_every_time = false;
		
		// if (max > MIN_PS_FOR_NBSEARCH && !isyielding){ //First time yielding, data has not been cleared from first search
			// ClearNbData(); 
			// MainNeighbourSearch/*_Ext*/();
			// isyielding  = true ;
		// }
		// if ( max > MIN_PS_FOR_NBSEARCH || isfirst || check_nb_every_time){	//TO MODIFY: CHANGE
			// if ( ts_i == 0 ){

				// if (m_isNbDataCleared){
					// MainNeighbourSearch/*_Ext*/();
          // //if (contact) SaveContNeighbourData();
					
	
				// }// ts_i == 0				
				
			// }
		
    // } //( max > MIN_PS_FOR_NBSEARCH || isfirst ){	//TO MODIFY: CHANGE

		// //NEW, gradient correction
			if (isfirst) {
				if (gradKernelCorr){
          CalcGradCorrMatrix();	}
				
			}		
		//std::cout << "neighbour_time (chrono, clock): " << clock_time_spent << ", " << neighbour_time.count()<<std::endl;
		
		GeneralBefore(*this);
		PrimaryComputeAcceleration();
    
       
    //CalcDensInc(); //TODO: USE SAME KERNEL?
    CalcRateTensorsDens();
    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (size_t i=0; i<Particles.Size(); i++){
      //Particles[i]->UpdateDensity_Leapfrog(deltat);
      Particles[i]->Density += dt*Particles[i]->dDensity;
    }    
    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (size_t i=0; i<Particles.Size(); i++){
      //Particles[i]->Mat2Leapfrog(deltat); //Uses density  
      Particles[i]->CalcStressStrain(deltat); //Uses density  
    }   

    CalcAccel(); //Nor density or neither strain rates    
    
    //BEFORE
    Vec3_t du,dv;
    
    #pragma omp parallel for schedule (static) private(du,dv) num_threads(Nproc)
    for (size_t i=0; i<Particles.Size(); i++){
      dv = Particles[i]->a*dt;
      du = Particles[i]->v*dt;
      Particles[i]->Displacement += du + dv*dt/2.;
      Particles[i]->x += du + dv*dt/2.;
    }      

    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (size_t i=0; i<Particles.Size(); i++){
      Particles[i]->v += Particles[i]->a*dt;
    }   
    GeneralAfter(*this);
    
		steps++;
		//cout << "steps: "<<steps<<", time "<< Time<<", tout"<<tout<<endl;
		// output
		if (Time>=tout){
			if (TheFileKey!=NULL) {
				String fn;
				fn.Printf    ("%s_%04d", TheFileKey, idx_out);
				WriteXDMF    (fn.CStr());
				//fn.Printf    ("%s_%.5f", TheFileKey, Time);
				WriteCSV    (fn.CStr());

			}
			idx_out++;
			tout += dtOut;
			std::cout << "\nOutput No. " << idx_out << " at " << Time << " has been generated" << std::endl;
			std::cout << "Current Time Step = " <<deltat<<std::endl;
			cout << "Max plastic strain: " <<max<< "in particle" << imax << endl;
			cout << "Max Displacements: "<<max_disp<<endl;
		}

		Time += deltat;

		
		// if (max>MIN_PS_FOR_NBSEARCH){	//TODO: CHANGE TO FIND NEIGHBOURS
			// if ( ts_i == (ts_nb_inc - 1) ){
				// ClearNbData();
			// }

			// ts_i ++;
			// if ( ts_i > (ts_nb_inc - 1) ) 
				// ts_i = 0;
		
		// }
    
    // if (Particles[0]->FirstStep)
    // for (size_t i=0; i<Particles.Size(); i++){
      // Particles[i]->FirstStep = false;
    // }
		if (isfirst) isfirst = false;
	
	}
	

	of.close();
	
	std::cout << "\n--------------Solving is finished---------------------------------------------------" << std::endl;

}

// THIS IS LIKE THE FRASER ALGORITHM
inline void Domain::SolveDiffUpdateModVerlet (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx) {
	std::cout << "\n--------------Solving---------------------------------------------------------------" << std::endl;

	size_t idx_out = 1;
	double tout = Time;

	//Initializing adaptive time step variables
	deltat = deltatint = deltatmin	= dt;
	

	InitialChecks();
	CellInitiate();
	ListGenerate();
	PrintInput(TheFileKey);
	TimestepCheck();
	WholeVelocity();
	

	//Initial model output
	if (TheFileKey!=NULL) {
		String fn;
		fn.Printf    ("%s_Initial", TheFileKey);
		WriteXDMF    (fn.CStr());
		std::cout << "\nInitial Condition has been generated\n" << std::endl;
	}
	

	unsigned long steps=0;
	unsigned int first_step;
	
	int ts_i=0;

	bool isfirst = true;
	bool isyielding = false;
  
  cout << "Nb search"<<endl;
  ClearNbData();
  cout << "cleared"<<endl;
  MainNeighbourSearch();
  //SaveNeighbourData();				//Necesary to calulate surface! Using Particle->Nb (count), could be included in search
  cout << "Done"<<endl;
  //CalculateSurface(1);				//After Nb search	
	//ClearNbData();
	
	//Print history
	std::ofstream of("History.csv", std::ios::out);
  of << "Displacement, pl_strain, eff_strain_rate, sigma_eq, sigmay, contforcesum"<<endl;
  
  bool check_nb_every_time = false;

  cout << "Main Loop"<<endl;
  int ct =30;
	while (Time<=tf && idx_out<=maxidx) {
    
		StartAcceleration(Gravity);

		double max = 0;
		int imax;
		#pragma omp parallel for schedule (static) num_threads(Nproc)	//LUCIANO//LIKE IN DOMAIN->MOVE
		for (int i=0; i<Particles.Size(); i++){
			if (Particles[i]->pl_strain > max){
        omp_set_lock(&dom_lock);
				max= Particles[i]->pl_strain;
        omp_unset_lock(&dom_lock);
				imax=i;
			}
		}

    Vec3_t max_disp = Vec3_t(0.,0.,0.);
		for (int i=0; i<Particles.Size(); i++){
      for (int j=0;j<3;j++)
        if (Particles[i]->Displacement[j]>max_disp[j]){
          max_disp[j] = Particles[i]->Displacement [j];
          imax=i;
			}
		}
    
    // ATTENTION! COULD BE LARGE DISPLACEMENTS AND SMALL STRAINS 
    //EXAMPLE COMPRESSION WITH NO FRICTION, SO CONTACTS NBs SHOULD BE RECALCULATED
    // if (norm(max_disp) > 0.1 * hmax){
      // if (!check_nb_every_time)
        // cout << "Checking Nb Every step now."<<endl;
      // check_nb_every_time = true;
    // }
    // else 
      // check_nb_every_time = false;
		
		// if (max > MIN_PS_FOR_NBSEARCH && !isyielding){ //First time yielding, data has not been cleared from first search
			// ClearNbData(); 
			// MainNeighbourSearch/*_Ext*/();
			// isyielding  = true ;
		// }
		// if ( max > MIN_PS_FOR_NBSEARCH || isfirst || check_nb_every_time){	//TO MODIFY: CHANGE
			// if ( ts_i == 0 ){

				// if (m_isNbDataCleared){
					// MainNeighbourSearch/*_Ext*/();
          // //if (contact) SaveContNeighbourData();
					
	
				// }// ts_i == 0				
				
			// }
		
    // } //( max > MIN_PS_FOR_NBSEARCH || isfirst ){	//TO MODIFY: CHANGE

		// //NEW, gradient correction
			if (isfirst) {
				if (gradKernelCorr){
          CalcGradCorrMatrix();	}
				
			}		
		//std::cout << "neighbour_time (chrono, clock): " << clock_time_spent << ", " << neighbour_time.count()<<std::endl;
		
		GeneralBefore(*this);
		PrimaryComputeAcceleration();
    
       
    //CalcDensInc(); //TODO: USE SAME KERNEL?
    CalcRateTensorsDens();
    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (size_t i=0; i<Particles.Size(); i++){
      //Particles[i]->UpdateDensity_Leapfrog(deltat);
      if (ct==30){
        Particles[i]->Densityb		= Particles[i]->Density;
        Particles[i]->Density			+=dt*Particles[i]->dDensity;
      } else {
        //Particles[i]->Density += dt*Particles[i]->dDensity;
        double dens	= Particles[i]->Density;
        Particles[i]->Density		= Particles[i]->Densityb + 2.0*dt*Particles[i]->dDensity;
        Particles[i]->Densityb	= dens;
      }
    }    
    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (size_t i=0; i<Particles.Size(); i++){
      //Particles[i]->Mat2Leapfrog(deltat); //Uses density  
      Particles[i]->CalcStressStrain(deltat); //Uses density  
    }   

    CalcAccel(); //Nor density or neither strain rates    
    
    //BEFORE
    Vec3_t du;
    #pragma omp parallel for schedule (static) private(du) num_threads(Nproc)
    for (size_t i=0; i<Particles.Size(); i++){
      du = (Particles[i]->v+Particles[i]->VXSPH) * dt + Particles[i]->a*dt*dt*0.5;
      Particles[i]->Displacement += du;
      Particles[i]->x += du;
    }      
    GeneralAfter(*this);
    
    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (size_t i=0; i<Particles.Size(); i++){
      if (ct == 30){
        Particles[i]->vb	= Particles[i]->v;
        Particles[i]->v	+=dt*Particles[i]->a;
      } else {
        Vec3_t temp;
        temp	= Particles[i]->v;
        Particles[i]->v		= Particles[i]->vb + 2*dt*Particles[i]->a;
        Particles[i]->vb		= temp;        
      }
    }   
    GeneralAfter(*this);
    


    
		steps++;
    if (ct==30) ct=0; else ct++;
		//cout << "steps: "<<steps<<", time "<< Time<<", tout"<<tout<<endl;
		// output
		if (Time>=tout){
			if (TheFileKey!=NULL) {
				String fn;
				fn.Printf    ("%s_%04d", TheFileKey, idx_out);
				WriteXDMF    (fn.CStr());
				//fn.Printf    ("%s_%.5f", TheFileKey, Time);
				WriteCSV    (fn.CStr());

			}
			idx_out++;
			tout += dtOut;
			std::cout << "\nOutput No. " << idx_out << " at " << Time << " has been generated" << std::endl;
			std::cout << "Current Time Step = " <<deltat<<std::endl;
			cout << "Max plastic strain: " <<max<< "in particle" << imax << endl;
			cout << "Max Displacements: "<<max_disp<<endl;
		}

		Time += deltat;

		
		// if (max>MIN_PS_FOR_NBSEARCH){	//TODO: CHANGE TO FIND NEIGHBOURS
			// if ( ts_i == (ts_nb_inc - 1) ){
				// ClearNbData();
			// }

			// ts_i ++;
			// if ( ts_i > (ts_nb_inc - 1) ) 
				// ts_i = 0;
		
		// }
    
    // if (Particles[0]->FirstStep)
    // for (size_t i=0; i<Particles.Size(); i++){
      // Particles[i]->FirstStep = false;
    // }
		if (isfirst) isfirst = false;
	
	}
	

	of.close();
	
	std::cout << "\n--------------Solving is finished---------------------------------------------------" << std::endl;

}


}; // namespace SPH
