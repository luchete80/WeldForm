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

#ifndef SPH_DOMAIN_H
#define SPH_DOMAIN_H

#include <stdio.h>    // for NULL
#include <algorithm>  // for min,max

#include <hdf5.h>
#include <hdf5_hl.h>

#include <omp.h>

#include "Particle.h"
#include "Functions.h"
#include "Boundary_Condition.h"

//#ifdef _WIN32 /* __unix__ is usually defined by compilers targeting Unix systems */
#include <sstream>
//#endif
#include <sstream>
#include <string>
#include <cmath>

#include "Mesh.h"
#include "Plane.h"

#include <fstream>

//#define NONLOCK_SUM 
#define MAX_NB_PER_PART 100

enum domain_bid_type {PlaneStress=0, PlaneStrain=1, AxiSymmetric =2, AxiSymm_3D =2};

//C++ Enum used for easiness of coding in the input files
enum Kernels_Type { Qubic_Spline=0, Quintic=1, Quintic_Spline=2 ,Hyperbolic_Spline=3};
enum Viscosity_Eq_Type { Morris=0, Shao=1, Incompressible_Full=2, Takeda=3 };
enum Gradient_Type { Squared_density=0, Multiplied_density=1 , Randles_Libersky = 2};

enum Damage_Model { None=0, Rankine=1 };

//TODO: Add from which it depends
enum Function_Type { Constant=0, Linear=1, Multilinear=2};

enum Friction_Type{Fr_Sta=0,Fr_Dyn,Fr_StaDyn,Fr_Bound};
enum Contact_Alg{Fraser=0, Wang, Seo, Zhan, LSDyna};

namespace SPH {
  
  //template <typename T>
  struct amplitude {
	int id;
	std::vector <double> time; ///set??
	//std::vector <T> value;
  //T getValAtTime(const double &t){
  std::vector <double> value;
  double getValAtTime(const double &t){
    double ret;
    //assumed ordered
    int i=time.size()-1;
    bool end = false;
    while (!end){
      i--;
      if (i==0 || t > time[i]) end =true;
    }
    ret =  value[i]+ (value[i+1] - value[i])/(time[i+1] - time[i])*(t-time[i]);
    return ret;
    
  }
	//std::map;
};

class NastranVolReader;


struct boundaryCondition {
	int 	zoneId;
	int 	type;	// ENUM TYPE Velocity, Force, Temperature
	bool 	free;	//is necessary??
	int 	valueType;		//0: Constant, 1 amplitude table
	Vec3_t value;       //If constant
  Vec3_t value_ang;       //Angular value
	int 	ampId;			//if valuetype == 1
	double 	ampFactor;		//if valuetype == 1
};

class Domain
{
public:
	typedef void (*PtVel) (Vec3_t & position, Vec3_t & Vel, double & Den, Boundary & bdry);
	typedef void (*PtOut) (Particle * Particles, double & Prop1, double & Prop2,  double & Prop3);
	typedef void (*PtDom) (Domain & dom);
    // Constructor
    Domain();

    // Destructor
    ~Domain();

    // Domain Part
    void AddSingleParticle	(int tag, Vec3_t const & x, double Mass, double Density, double h, bool Fixed);		//Add one particle
    void AddBoxLength				(int tag, Vec3_t const &V, double Lx, double Ly, double Lz,double r, double Density,
																	double h,int type, int rotation, bool random, bool Fixed);									//Add a cube of particles with a defined dimensions

	void AddCylinderLength(int tag, Vec3_t const & V, double Rxy, double Lz, 
									double r, double Density, double h, bool Fixed, bool ghost = false);

  void AddCylUniformLength(int tag, double Rxy, double Lz, 
																				double r, double Density, double h, double ang = 2 * M_PI, int rows=1);
                                        
	//Cylinder Slice 
	void AddXYSymCylinderLength(int tag, double Rxy, double Lz, 
								double r, double Density, double h, bool Fixed, bool symlength = false);
  void AddCylSliceLength(int tag, double alpha, double Rxy, double Lz, 
																				double r, double Density, double h);
  // void AddCylUniformLength    (int tag, double Rxy, double Lz, 
																				// double r, double Density, double h);                                  
	void AddTractionProbeLength(int tag, Vec3_t const & V, double Rxy, double Lz_side,
											double Lz_neckmin,double Lz_necktot,double Rxy_center,
											double r, double Density, double h, bool Fixed);
											
	void Calculate3DMass(double Density);
	void Add3DCubicBoxParticles(int tag, Vec3_t const & V, double Lx, double Ly, double Lz, 
									double r, double Density, double h);


    void AddBoxNo						(int tag, Vec3_t const &V, size_t nx, size_t ny, size_t nz,double r, double Density,
																	double h,int type, int rotation, bool random, bool Fixed);									//Add a cube of particles with a defined numbers
    void DelParticles				(int const & Tags);					//Delete particles by tag
    void CheckParticleLeave	();													//Check if any particles leave the domain, they will be deleted

    void YZPlaneCellsNeighbourSearch(int q1);						//Create pairs of particles in cells of XZ plan
    void MainNeighbourSearch				();									//Create pairs of particles in the whole domain
    void MainNeighbourSearch_CNS    (const double &r);  //NEW; ALLOWS TO SAVE DATA BY PARTICLE NBS (AND NOT LOCKING DOMAIN)
		void MainNeighbourSearch_Ext		();									//Create pairs of particles in the whole domain
		int AvgNeighbourCount						();									//Create pairs of particles in the whole domain
    
    void InitReductionArraysOnce();
    inline void ResetReductionArrays();
    inline void CalcPairPosList();                             //Calculate position list for every particle ipl/jpl[NProc][particle]
    inline void CalcRefTable();
    //For new reduction method
    inline void AccelReduction();
    inline void RateTensorsReduction();
    inline void DensReduction();
    void CheckParticlePairs(const int &i);
		
		void SaveNeighbourData();
		void SaveContNeighbourData();
    
    void ApplyAxiSymmBC(int bc_1 = -1, int bc_2 = -1); //Apply to all particles or only to BCs. If Not all (!=-1), 
		
    void StartAcceleration					(Vec3_t const & a = Vec3_t(0.0,0.0,0.0));	//Add a fixed acceleration such as the Gravity
    void PrimaryComputeAcceleration	();									//Compute the solid boundary properties
    void LastComputeAcceleration		();									//Compute the acceleration due to the other particles
    void CalcForce2233	(Particle * P1, Particle * P2);		//Calculates the contact force between soil-soil/solid-solid particles
    void CalcAccel();		//NEW, ONLY CALCULATES ACCELERATION; IN ORDER TO ALTERNATE AND NOT CALCULATE Density at same place
    
    ///// TEST FUNCTIONS FOR PARTICLE PARALLELIZATION
    inline void CalcAccelPP(); //ONLY FOR TESTING, PARALLELIZATION BY PARTICLE
    inline void CalcAccelPair(Particle * P1, Particle * P2);
    inline void CalcRateTensorsPair (Particle *P1, Particle *P2);
    inline void CalcTensorsPP();
    inline void CalcDensPP();
    inline void CalcDensIncPairs(Particle *P1, Particle *P2);
    
    void CalcDensInc();
    void CalcRateTensors();
    void CalcForceSOA(int &i,int &j) ;
    void Move						(double dt);										//Move particles

  void Solve					(double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx);		///< The solving function
  void SolveDiffUpdateLeapfrog(double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx);
  void SolveDiffUpdateKickDrift (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx);
  void SolveDiffUpdateLeapFrog (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx);
  void SolveDiffUpdateVerlet  (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx);
  void SolveDiffUpdateModEuler (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx);
  void SolveDiffUpdateFraser (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx);
  void Solve_orig 			(double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx); 
  void Solve_orig_Ext 	(double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx); 
	
  // THERMAL //
  void ThermalSolve			(double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx);		///< The solving function
	void ThermalStructSolve (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx); //Coupled Thermal Structural
	void ThermalSolve_wo_init	(double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx);		///< The solving function
  inline void ThermalCalcs(const double &dt);
  
  //inline void CalcContactFrictionHeat();
  
  
  void AddFixedMassScaling (const double &factor);
  inline void UpdateSmoothingLength();
  //inline void UpdateSmoothingLength_Pairs();
	
	inline void CalcThermalExpStrainRate();
	inline void CalcPlasticWorkHeat(const double &dt);

	inline void CalcKinEnergyEqn();
  inline void CalcIntEnergyEqn();

    void Solve_wo_init (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx);		///< The solving function	
	//void Step(double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx);
	
    void CellInitiate		();															//Find the size of the domain as a cube, make cells and HOCs
    void ListGenerate		();															//Generate linked-list
    void CellReset			();															//Reset HOCs and particles' LL to initial value of -1
	
	void ClearNbData();	
	
    void WriteXDMF			(char const * FileKey);					//Save a XDMF file for the visualization
    void WriteCSV				(char const * FileKey);					//Save a XDMF file for the visualization
    
    void ReadXDMF			(char const * FileKey);	        //NEW, FOR RESTART

    void InFlowBCLeave	();
    void InFlowBCFresh	();
    void WholeVelocity	();

	void Kernel_Set									(Kernels_Type const & KT);
	void Viscosity_Eq_Set						(Viscosity_Eq_Type const & VQ);
	void Gradient_Approach_Set			(Gradient_Type const & GT);
	
	//Thermal Solver
	void CalcTempInc 		(); 		//LUCIANO: Temperature increment
	void CalcTempIncSOA (); 		//LUCIANO: Temperature increment
  inline double CalcTempIncPair (Particle *P1, Particle *P2);  //PER PARTICLE; ONLY FOR TEST
  inline void CalcTempIncPP();
	inline void CalcConvHeat ();
	inline void CalcConvHeatSOA();
	inline void CalcGradCorrMatrix();	//BONET GRADIENT CORRECTION
	inline void CalcGradCorrMixedMatrix();	//BONET GRADIENT CORRECTION
  
  inline void CorrectVelAcc();  //For axil symm particles
	
	inline void MoveGhost();
  inline void PropGhost();
	const double & getStepSize()const {return deltat;};

	
	
	/////////////////////// CONTACT /////////////////////////////
	void AddTrimeshParticles(TriMesh *mesh, const float &hfac, const int &id);
  void CalculateSurface(const int &id = 1);
  
  void GenerateSPHMesh(const int &tag, NastranVolReader &nr,double Density, double hfac);
	
  inline void CalcContactForces();
  inline void CalcContactForcesLS();
  inline void CalcContactForcesAnalytic();
  inline void CalcContactForces2(); //Position criteria, SEO Contact detection
  inline void CalcContactForcesWang();
  inline void CalcContactInitialGap();
  inline void UpdateContactParticles();  //Update position, velocity and normals FROM MESH
  
  inline void UpdateFrictionCoeff();
  
  //////////////////////////// ENERGY
  double kin_energy_sum, int_energy_sum;
  double mass_scaling_factor;
  
  bool contact_mesh_auto_update;
  inline void ContactNbSearch();	//Performed AFTER neighbour search
	std::vector<int> contact_surf_id;						//particles id from surface

	double contact_force_factor;
	double friction;
  double friction_sta, friction_dyn;
  double friction_m, friction_b;  //linear coeff TODO: MOVE TO FUNCTION
  
  Function_Type friction_function;
  // double friction_sta;
  // double friction_kin;
  Friction_Type fric_type;
  
  //Profiler things
  double m_contact_forces_time;
  double m_clock_begin;
  double m_forces_artifvisc_time;
  double m_forces_momentum_time;
  double m_forces_tensors_time;
  double m_forces_update_time;
	double contact_friction_work, plastic_work, ext_f_work;
  double pl_work_heat_frac; //Fration of plastic work converted to heat
  double ext_forces_work,ext_forces_work_step;
  
  int ts_nb_inc;
  
  int solid_part_count;
  //TEST
  //Forces calculation time spent
  double forces_tensile_inst_calctime, forces_stressstrain_calctime,forces_acc_calctime,
          forces_artif_visc_calctime;
	
	/////////////// MEMBERS //
    // Data
	//    std::vector< *Particle >				Particles; 	///< Array of particles
	Array <Particle*>				Particles; 	///< Array of particles
    double					R;		///< Particle Radius in addrandombox

		double					sqrt_h_a;				//Coefficient for determining Time Step based on acceleration (can be defined by user)
		double 					min_force_ts;		//min time step size due to contact forces
    int             manual_min_ts;
    int             min_ts_acc_part_id;
		
    int 					  Dimension;    	///< Dimension of the problem
    domain_bid_type  dom_bid_type;    


    double					MuMax;		///< Max Dynamic viscosity for calculating the timestep
    double					CsMax;		///< Max speed of sound for calculating the timestep
	double 					Vol;		///LUCIANO
	
		bool cont_heat_gen;     //Contact heat generation
    bool cont_heat_cond;    //Contact Heat conduction
    
    double contact_hc;      //Conductance coeff
    
    std::vector <double>    tot_cont_heat_cond;   //Total contact heat conductance
    double                  accum_cont_heat_cond; 
		std::vector <int> 			first_fem_particle_idx;			//The rest are ridig bodies
    int                     meshcount;
		
    /*Array<*/int/*>*/ 			id_free_surf;								//TODO: 
    Vec3_t					        Gravity;       	///< Gravity acceleration
    
    bool                    h_update;

    Vec3_t                 			TRPR;		///< Top right-hand point at rear of the domain as a cube
    Vec3_t                  			BLPF;           ///< Bottom left-hand point at front of the domain as a cube
    Vec3_t                  			CellSize;      	///< Calculated cell size according to (cell size >= 2h)
    int		                		CellNo[3];      ///< No. of cells for linked list
    double 					hmax;		///< Max of h for the cell size  determination
    Vec3_t                 			DomSize;	///< Each component of the vector is the domain size in that direction if periodic boundary condition is defined in that direction as well
    double					rhomax;

    int						*** HOC;	///< Array of "Head of Chain" for each cell

    bool					FSI;						///< Selecting variable to choose Fluid-Structure Interaction
		int						contact_type;		//0: no contact 1: node to surface 2: node 2 node
		bool					thermal_solver;
    Contact_Alg   contact_alg;
    
   std::vector <boundaryCondition> bConds;  //NEW, For BCond
	
	// BONET KERNEL CORRECTION
	bool 					gradKernelCorr;	
	
    double 					XSPH;		///< Velocity correction factor
    double 					InitialDist;	///< Initial distance of particles for Inflow BC

    double					AvgVelocity;	///< Average velocity of the last two column for x periodic constant velocity
	double 					getCellfac(){return Cellfac;}
  
  bool  model_damage; //Include or not a damage model, ENUM

	#ifdef __GNUC__
    size_t					Nproc;		///< No of threads which are going to use in parallel calculation
	#else
	int						Nproc;
	#endif
	omp_lock_t 					dom_lock;	///< Open MP lock to lock Interactions array
    Boundary					BC;
    PtOut					UserOutput;
    PtVel 					InCon;
    PtVel 					OutCon;
    PtVel 					AllCon;
    Vec3_t					DomMax;
    Vec3_t					DomMin;
    PtDom					GeneralBefore;	///< Pointer to a function: to modify particles properties before CalcForce function
    PtDom					GeneralAfter;	///< Pointer to a function: to modify particles properties after CalcForce function
    size_t					Scheme;		///< Integration scheme: 0 = Modified Verlet, 1 = Leapfrog
    
    double CFL;               ///FOR VELOCITY
    bool  nonlock_sum;        //Reduction by table instead of thread locking (Nishimura)

    Array<Array<std::pair<size_t,size_t> > >	SMPairs;
    Array<Array<std::pair<size_t,size_t> > >	NSMPairs;
    Array<Array<std::pair<size_t,size_t> > >	FSMPairs;
    
    ////// IF DAMAGE!!
    Array<Array<double>>                      dam_D;            //DAMAGE FACTOR!!(processor, pair)
    Array<Array<bool>>                        dam_pair;         //DAMAGED PAIR !!(processor, pair), USED TO CALC SURFACE
    
    ////// THERE ARE TWO INITIAL DISTANCES, CRACK
    Array<Array<double>>                      dam_r0;          //Initial distance BETWEEN PARTICLES; IF UNIFORM THIS ARRAY IS NOT ALLOCATED
    Array<Array<double>>                      dam_rf0;         //Initial failure distance
    double                                    dam_r0_unif;      //IF MESH IS UNIFORM rij_0 is also uniform
    bool                                      is_grid_uniform;  
    
		Array<Array<std::pair<size_t,size_t> > >	RIGPairs;	//Previous instance to contact pairs, because outer surface should be located
																												//based on original neighbours
		
    Array<Array<std::pair<size_t,size_t> > >	ContPairs;//Asuming same material
    
    //NEW: For parallel sum/reduction
    std::vector<size_t> first_pair_perproc;                   // Almost like pair count        
    int pair_count;                                                   //var names as stated as Nishimura (2015) ipl is njgi 
    std::vector < std::pair<int,int> >     pair_test,pair_ord;                    //OLY FOR TESTING
    //Array< Array <size_t> >               ilist_SM,jlist_SM;          // Size [Pairs] i and j particles of pair list [l], already flattened
    std::vector < size_t >                ipair_SM,jpair_SM;          //[Particles]// This is nb count for each particle i<j and j>i (called njgi) FLATTENED
    std::vector < size_t >                ipl_SM;                     // [Particles] position of link/pair (nb sum), called s_jgi in 1991 paper
    std::vector < std::vector <size_t>  > Aref;                      // Entry[Particle ,nb], indicates link
    std::vector < std::vector <size_t>  > Anei;                      //[Particles][MAX_NB_PER_PART] neighbiour list for j > i
    std::vector <Vec3_t>                  pair_force;
    std::vector <double>                  temp_force;
    std::vector <double>                  pair_densinc;
    std::vector <Mat3_t>                  pair_StrainRate;
    std::vector <Mat3_t>                  pair_RotRate;    
    
    Array< size_t > 				FixedParticles;
    Array< size_t >				FreeFSIParticles;
	double 	& getTime (){return Time;}		//LUCIANO
		bool 									update_contact_surface;

    Array<std::pair<size_t,size_t> >		Initial;
    Array<double>                       dam_initial;  //if damage
    Array<bool>                         dam_pair_initial;  //if damage
    
    Mat3_t I;
    String					OutputName[3];
	double T_inf;			//LUCIANO: IN CASE OF ONLY ONE CONVECTION TEMPERAURE
	
	iKernel m_kernel;
	bool					m_isNbDataCleared;
	bool						auto_ts;				//LUCIANO: Auto Time Stepping: VEL CRITERIA
	bool            auto_ts_acc;    
  bool            auto_ts_cont; 
  double          m_maxT, m_minT; //from step, tou output
  
  std::vector <SPH::Plane*> planes;
  std::vector <TriMesh*> trimesh; //ORIGINALLY
	//CONTACT 
	double PFAC, DFAC;		// Penalty and damping factors
	bool 		contact;
	double max_contact_force;
  
  double contact_force_sum;
  double contact_reaction_sum;  //Acceleration of contact particles previous to apply contact_force
  
  double m_scalar_prop;  //User Defined Domain Property
	Vec3_t max_disp;
    
  //ATTENTION: REDUNDANT, ghost pairs and reference
	Array<std::pair<size_t,size_t> > GhostPairs;	//If used
	
	/////////////////////// SOA (Since v0.4) ///////////////////////////////////
	Vec3_t **m_x,*m_v,*m_a;
	double **m_h;
	double **m_T, **m_Tinf, **m_kT, **m_hcT, **m_cpT, **m_dTdt;
	double **m_qconvT,**m_qT;	//thermal source terms 
	double **m_rho, **m_mass;
  //Mechanics
  double *sigma;
	double *strrate,*rotrate;//all flattened, six component (rotation rate is upped matrix component, since it is antisymm)
	double *shearstress,*shearstressa,*shearstressb;
	double *strain,*straina,*strainb;
	
  //////////////////////// NEW: IMPLICIT SOLVER FOR QUASI STATIC 
  inline void InitImplicitSolver();
  inline void CalcBMat();
  inline void CalcStiffMat();
  inline void CheckMinTSVel();
  inline void CheckMinTSAccel();
  inline void CalcDamage();
  
  int AssignZone(Vec3_t &start, Vec3_t &end, int &id);
	
  std::vector <SPH::amplitude> amps; ////maybe move to domain
  string filename;
  std::ofstream out_file;
  
  private:
		bool  Domain::CheckRadius(Particle* P1, Particle *P2);
		void Periodic_X_Correction	(Vec3_t & x, double const & h, Particle * P1, Particle * P2);		//Corrects xij for the periodic boundary condition
		void AdaptiveTimeStep				();		//Uses the minimum time step to smoothly vary the time step
    

		void PrintInput			(char const * FileKey);		//Print out some initial parameters as a file
		void InitialChecks	();		//Checks some parameter before proceeding to the solution
		void TimestepCheck	();		//Checks the user time step with CFL approach

		size_t					VisEq;					//Choose viscosity Eq based on different SPH discretisation
		size_t					KernelType;			//Choose a kernel
		size_t					GradientType;		//Choose a Gradient approach 1/Rho i^2 + 1/Rho j^2 or 1/(Rho i * Rho j)
		double 					Cellfac;				//Define the compact support of a kernel

		double					Time;    				//Current time of simulation at each solving step
		double					deltat;					//Time Step
    double					deltatmin;			//Minimum Time Step
    double					deltatint;			//Initial Time Step
		
		int 						cont_pairs;
		bool enable_th_exp;
		bool enable_plastic_heat_gen;
		void AllocateNbPair(const int &temp1, const int &temp2, const int &T);
    


		

};


}; // namespace SPH

#include "Interaction.cpp"
#include "InteractionTest.cpp"
#include "Domain.cpp"

#include "Neighbour.cpp"
//#include "Input.cpp"
#include "Output.cpp"
#include "InOutFlow.cpp"

#include "Thermal.cpp"
#include "ThermalStuct.cpp"
#include "Contact.cpp"

#include "ImplicitSolver.cpp"

#include "Damage.cpp"

#endif // SPH_DOMAIN_H
