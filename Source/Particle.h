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

#ifndef SPH_PARTICLE_H
#define SPH_PARTICLE_H

#include "matvec.h"
#include "Functions.h"

#include "Material.h" //Johnson Cook

#define TH_BC_NONE			0
#define TH_BC_CONVECTION	1
#define TH_BC_CONDUCTION	2

//TODO: PASS TO MATERIAL CLASS
#define BILINEAR				     0
#define HOLLOMON				1 //POWER LAW
#define JOHNSON_COOK		2
#define _GMT_	                           3

#include "Matrix.h"  /////ONLY FOR IMPLICIT SOLVER

//#define 	FLAT_TENSORS

#include "Plane.h" //ONLY FOR GHOST

enum Ghost_Type {Symmetric = 0, Periodic = 1, Mirror_XYZ = 2 };

class Plane;    //For analytical contact 
namespace SPH {

	class Particle
	{
	public:
		// Shepard density correction
		bool   	Shepard;	///< Shepard Filter for the density
		size_t	ShepardCounter;	///< Count number of contributing particles
		size_t	ShepardStep;	///< Cycle number for shepard correction
		double	ZWab;		///< Summation of mb/db*Wab for neighbour particles of the particle a (for Shepard filter)
		double	SumDen;		///< Summation of mb*Wab for neighbour particles of the particle a (for Shepard filter)

		bool   	IsFree;		///< Check the particle if it is free to move or not
		size_t	InOut;		///< Check the particle if it is in-flow or out-flow or not
		bool   	IsSat;		///< Check the particle if it is Saturated or not
		bool   	SatCheck;	///< Check the particle Saturation at each time step
		bool   	NoSlip;		///< No-Slip BC
		int			Material_model;	//

		int    	ID;		///< an Integer value to identify the particle set
		int 		ID_orig;
		int 	Thermal_BC;
		int    	Material;	///< an Integer value to identify the particle material type: 1 = Fluid, 2 = Solid, 3 = Soil
		Material_*	mat;	//NOT TO CONFUSE WITH MATERIAL NUMBER (TO BE DELETED)	
	
		Vec3_t	x;		///< Position of the particle n
		Vec3_t	vb;		///< Velocity of the particle n-1 (Modified Verlet)
		Vec3_t	va;		///< Velocity of the particle n+1/2 (Leapfrog)
		Vec3_t	v;		///< Velocity of the particle n+1
    Vec3_t	vc;   ///Contact velocity, for bonded contact
		Vec3_t	NSv;		///< Velocity of the fixed particle for no-slip BC
		Vec3_t	VXSPH;		///< Mean Velocity of neighbor particles for updating the particle position (XSPH)
    Vec3_t  v_max;
		Vec3_t	a;		///< Acceleration of the particle n
		bool		update;	/// UPDATE KERNELS AND NEIGHBOURS (IF DEP IS INCREMENTING IN LAST STEP IS SET TO TRUE)
		double  hfac;      //For smoothing length particle update

		size_t	PresEq;		///< Selecting variable to choose an equation of state
		double	Cs;		///< Speed of sound
		double	P0;		///< background pressure for equation of state
		double 	Pressure;	///< Pressure of the particle n+1
    
    Mat3_t  Bmat;   //////Randles & Libersky, CAMPBELL (2000), boundary correction matrix
    
    //Boundary corrections
    bool correct_vel_acc;
    bool is_fixed;
    bool is_boundary; // TO APPLY RANDLES & LIBERSKY COND
    
    double ps_energy; //PER UNIT VOLUME

		double	Density, Density_real;	///< Density of the particle n+1, Density_real is the real density used in Axisymm calcs
		double 	Densitya;	///< Density of the particle n+1/2 (Leapfrog)
		double 	Densityb;	///< Density of the particle n-1 (Modified Verlet)
		double 	dDensity;	///< Rate of density change in time based on state equations n
    double  etaDens; //Desity in axi-symm problems
		double 	RefDensity;	///< Reference Density of Particle
		double 	FPMassC;	///< Mass coefficient for fixed particles to avoid leaving particles
		double 	Mass;		///< Mass of the particle
		Vec3_t	Displacement;	///< Density of the particle n+1
    Vec3_t  x_prev;           //ONLY FOR SEO CONTACT (contact in current, not predicted, step)z
    double  friction_hfl; //Surface 
    double  cshearabs;   // Contact shear stress module, for comparison
    
		Mat3_t	StrainRate;	///< Global shear Strain rate tensor n
		Mat3_t	RotationRate;	///< Global rotation tensor n
		double	ShearRate;	///< Global shear rate for fluids
		double	SBar;		///< shear component for LES
    double  eff_strain_rate;
                
		double strrate[6], rotrate[3];	//Null Diagonal for rot rate
		Mat3_t	ShearStress;	///< Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor) n+1
		Mat3_t	ShearStressa;	///< Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor) n+1/2 (Leapfrog)
		Mat3_t	ShearStressb;	///< Deviatoric shear stress tensor (deviatoric part of the Cauchy stress tensor) n-1 (Modified Verlet)

		Mat3_t	Sigma;		///< Cauchy stress tensor (Total Stress) n+1

		Mat3_t	Sigmaa;		///< Cauchy stress tensor (Total Stress) n+1/2 (Leapfrog)
		Mat3_t	Sigmab;		///< Cauchy stress tensor (Total Stress) n-1 (Modified Verlet)
		
		double Sigma_eq;	//Von Mises
		
		////////////////// PLASTIC THINGS
		Mat3_t	Strain;							///< Total Strain n+1
		Mat3_t	Straina;					///< Total Strain n+1/2 (Leapfrog)
		Mat3_t	Strainb;					///< Total Strain n-1 (Modified Verlet)
		Mat3_t  Strain_pl;				//// Plastic Strain
		Mat3_t  Strain_pl_incr;		//// Plastic Strain - INTERNAL, JUST FOR PLASTIC THERMAL HEAT GEN CALCULATION
    
    Matrix  m_B;              //B matrix for strain (IMPLICIT SOLVER)
		
    Plane  *plane;
		
		double 	pl_strain,delta_pl_strain;	//Accum and incremental Effective (Von Mises) plastic strain 
		
		// BONET GRADIENT CORRECTION MATRIX
		Mat3_t	gradCorrM;

		Mat3_t	TIR;		///< Tensile Instability stress tensor R
		double	TI;		///< Tensile instability factor
		double	TIn;		///< Tensile instability power
		double 	TIInitDist;	///< Initial distance of particles for calculation of tensile instability

		double 	Alpha;		///< Dynamic viscosity coefficient of the fluid particle
		double 	Beta;		///< Dynamic viscosity coefficient of the fluid particle
		double 	Mu;		///< Dynamic viscosity coefficient of the fluid particle
		double 	MuRef;		///< Reference Dynamic viscosity coefficient
		double 	T0;		///< Yield stress for Bingham fluids
		double 	m;		///< Normalization value for Bingham fluids
		size_t	VisM;		///< Non-Newtonian viscosity method
		bool		LES;		///< Large eddy simulation using sub-particle scale
		double	CSmag;		///< Coefficient of Smagorinsky-Lilly model

		double 	G;		///< Shear modulus //TODO: Move To Material
		double 	K;		///< Bulk modulus  //TODO: Move To Material
		double	Sigmay;		///< Tensile yield stress
		size_t	Fail;		///< Failure criteria
		
		double 	Ep, Et;		//TODO: Move To Material
		double 	Et_m;	//This not make sens to be a material member since is instantaneous 

		double	V;		///< Volume of a particle

		double 	h,hmin,hmax,hini;		///< Smoothing length of the particle
		int    	LL;		///< Linked-List variable to show the next particle in the list of a cell
		int    	CC[3];		///< Current cell No for the particle (linked-list)
		int		ct;		///< Correction step for the Modified Verlet Algorithm
		double	SumKernel;	///< Summation of the kernel value for neighbour particles
		bool	FirstStep;	///< to initialize the integration scheme
		bool	ThermalFirstStep;	///< to initialize the integration scheme
		bool 	print_history; 		//TODO: CREATE AN ENUM TO VARIABLES
    
    Ghost_Type ghost_type;
		
		
		//LUCIANO: THERMAL PROPERTIES
		double T,k_T,cp_T,dTdt;			// Temperature, avoid permeability	
		double Ta,Tb;						//Temperature (t-1) for leapfrog
		double q_source;
		double q_conv,T_inf,h_conv;				//Different heat source terms
    double q_cont_conv;
    
    double mcp_t;                   //mass * cp, ONLY FOR CONTACT SOLID PARTICLES 
    
		double q_plheat;				//Plastic Work Heat generation
		double th_exp;		//Constant
		
		double dkin_energy_dt;		//
		double dint_energy_dt;		//
		double kin_energy;
		double int_energy;
		
		double q_fric_work;
		
		int 	Nb;
		int 	ContNb;
		
		omp_lock_t my_lock;		///< Open MP lock
    
    bool not_write_surf_ID;
		
		////////////////// CONTACT //////////////////////
		double 	cont_stiff;
    double  delta_cont; // penetration
		Vec3_t	contforce;
    Vec3_t	tgforce;
    int     mesh;       //Corresponding rigid mesh
		
		Vec3_t	tgdir;			//TEMP: FO DEBUG
		/////////////////////// SURFACE NORMALS /////////////////////
		Vec3_t normal;	
		bool is_surface;		// BUT A SPECIFIC ID WILL BE ADDED
		int element;				// Element index (if comes from FEM)
		
    //SYMMETRY AND GHOST
    int     inner_mirr_part;       //NEW: SYMMETRIC PARTICLE 
    bool    is_ghost;      
    int     ghost_plane_axis;      //Assuming is cartesian
		bool    impose_vel;
    
    SPH::Plane *plane_ghost;
    
    //// DAMAGE (NOT USED BY RANKINE METHOD)
    double dam_D;
			
		// Constructor
		Particle						(int Tag, Vec3_t const & x0, Vec3_t const & v0, double Mass0, double Density0, double h0, bool Fixed=false);

		// Methods
		void Move						(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin,size_t Scheme, Mat3_t I);	///< Update the important quantities of a particle
		void Move_MVerlet		(Mat3_t I, double dt);										///< Update the important quantities of a particle
		void Move_Verlet		(Mat3_t I, double dt);		//LUCIANO
		void Move_Leapfrog	(Mat3_t I, double dt);										///< Update the important quantities of a particle
    //These two are set in order to alternate update (Randles Libersky 1996)
    void UpdateDensity_Leapfrog(double dt);
    void UpdateVelPos_Leapfrog(double dt);
    
		void translate			(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin);
		void Mat2Verlet			(double dt);
		void Mat2MVerlet		(double dt);
		void TempCalcLeapfrog	(double dt);
		void Mat2Leapfrog		(double dt);
    void CalcStressStrain (double dt);
		void PlasticHeatTest	();
		void CalcPlasticWorkHeat(const double &dt);
		void CalcThermalExpStrainRate();
		void CalculateEquivalentStress();
		void Move_Euler (Mat3_t I, double dt);
		void Mat2Euler(double dt);
		void CalcIntEnergyEqn();
    inline void LimitVel();
    
    //Implicit Solver
    inline void AddBMat(const Vec3_t &);
		
		bool is_axisymm;  //FOR EOS CALC, should alter density

	};
}; // namespace SPH

#include "Particle.cpp"

#endif //SPH_PARTICLE_H
