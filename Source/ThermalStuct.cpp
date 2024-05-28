#include "Domain.h"
#include <vector>
using namespace std;

namespace SPH {

inline void Domain::CalcThermalExpStrainRate(){
	#pragma omp parallel for schedule (static) num_threads(Nproc)
		for (int p=0;p<Particles.Size();p++){
			Particles[p]->CalcThermalExpStrainRate();	//Add Thermal expansion Strain Rate Term
		}		
}

inline void Domain::CalcPlasticWorkHeat(const double &dt){
  if (pl_heating){
    double f;
    double plw, pl_sum = 0.;
    #pragma omp parallel for schedule (static) num_threads(Nproc) private (f)
    for (int p=0;p<Particles.Size();p++){
      if (p < solid_part_count){
        f = Particles[p]->Mass / Particles[p]->Density; //Vol
        Particles[p]->CalcPlasticWorkHeat(dt);	//Add Thermal expansion Strain Rate Term
        
      
        // plw = f * Particles[p]->q_plheat;
        // //#pragma omp atomic
        // omp_set_lock(&dom_lock);            
          // pl_sum += plw;     
        // omp_unset_lock(&dom_lock);	
      }
    }
    //Serialized
    for (int p=0;p<solid_part_count;p++)
        pl_sum += Particles[p]->Mass / Particles[p]->Density *Particles[p]->q_plheat;
      
    plastic_work += pl_sum * dt;  
    //cout << "Max plastic heat gen: "<<max<<endl;
  }
}

inline void Domain::ThermalStructSolve (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx) {
	std::cout << "\n--------------Solving---------------------------------------------------------------" << std::endl;

	size_t idx_out = 1;
	double tout = Time;

	//Initializing adaptive time step variables
	deltat = deltatint = deltatmin	= dt;
	
	auto start_whole = std::chrono::steady_clock::now();

	InitialChecks();
	CellInitiate();
	ListGenerate();
	PrintInput(TheFileKey);
	TimestepCheck();
	WholeVelocity();
	
	std::chrono::duration<double> total_time,neighbour_time;
	
	clock_t clock_beg;
	
	double neigbour_time_spent_per_interval=0.;
	
	double clock_time_spent,pr_acc_time_spent,acc_time_spent;
	clock_time_spent=acc_time_spent=0.;
	
	//Initial model output
	if (TheFileKey!=NULL) {
		String fn;
		fn.Printf    ("%s_Initial", TheFileKey);
		WriteXDMF    (fn.CStr());
		std::cout << "\nInitial Condition has been generated\n" << std::endl;
	}

	int ts_nb_inc=3;	// Always > 0
	int ts_i=0;

	bool isfirst = true;
	bool isyielding = false;

	//In case of contact this must be SURFACE particles
	//TODO, REMOVE so many nb search
	if (contact){
		MainNeighbourSearch();
		SaveNeighbourData();				//Necesary to calulate surface! Using Particle->Nb (count), could be included in search
		CalculateSurface(1);				//After Nb search			
	}
	
	ClearNbData();

	while (Time<tf && idx_out<=maxidx) {

		StartAcceleration(Gravity);
		// //if (BC.InOutFlow>0) InFlowBCFresh();
		// auto start_task = std::chrono::system_clock::now();
		

		double max = 0;
		int imax;
		#pragma omp parallel for schedule (static) num_threads(Nproc)	//LUCIANO//LIKE IN DOMAIN->MOVE
		for (int i=0; i<Particles.Size(); i++){
			if (Particles[i]->pl_strain>max){
				max= Particles[i]->pl_strain;
				imax=i;
			}
		}
		
		if (max > MIN_PS_FOR_NBSEARCH && !isyielding){ //First time yielding, data has not been cleared from first search
			ClearNbData();

			MainNeighbourSearch();
			isyielding  = true ;
		}
		if ( max > MIN_PS_FOR_NBSEARCH || isfirst ){	//TO MODIFY: CHANGE
			if ( ts_i == 0 ){
				clock_beg = clock();
				if (m_isNbDataCleared)
					//cout << "Performing Nb search"<<endl;
					MainNeighbourSearch();

				neigbour_time_spent_per_interval += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
			}
			isfirst = false;
		}
		
		auto start_task = std::chrono::system_clock::now();
		clock_beg = clock();
		clock_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		auto end_task = std::chrono::system_clock::now();
		 neighbour_time = /*std::chrono::duration_cast<std::chrono::seconds>*/ (end_task- start_task);

		GeneralBefore(*this);
		clock_beg = clock();
		PrimaryComputeAcceleration();
		pr_acc_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		clock_beg = clock();
		LastComputeAcceleration();	//Define Strain Rate, and accelertion from prev sigma
		acc_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;

		CalcConvHeat();
		CalcPlasticWorkHeat(deltat);
		CalcTempInc();
		CalcThermalExpStrainRate();	//Add Thermal expansion Strain Rate Term		
		
		GeneralAfter(*this);

		if (auto_ts)
			AdaptiveTimeStep();
		
			// if (Time>0.1)
				// for (size_t i=0; i<Particles.Size(); i++){
				// cout <<"StrainRate"<<Particles[i]->StrainRate<<endl;}
		Move(deltat);		 //Calculates Sigma, u && v from accel
			//if (Time>0.1)
				// for (size_t i=0; i<Particles.Size(); i++){
				// cout <<"ShearStress"<<Particles[i]->ShearStress<<endl;}
		Time += deltat;
		
		//std::cout << "neighbour_time (chrono, clock): " << clock_time_spent << ", " << neighbour_time.count()<<std::endl;
		clock_beg = clock();

		//TODO Add 
		double Tmin=1000.;
		double Tmax=0.0;
		for (size_t i=0; i<Particles.Size(); i++){
			Particles[i]->TempCalcLeapfrog(deltat);
			if (Particles[i]->T > Tmax)
				Tmax=Particles[i]->T;
			if (Particles[i]->T < Tmin)
				Tmin=Particles[i]->T;
		}
		//std::cout << "Max temp: "<< max << std::endl;

			
		acc_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;

		// output
		if (Time>=tout){
			if (TheFileKey!=NULL) {
				String fn;
				fn.Printf    ("%s_%04d", TheFileKey, idx_out);
				WriteXDMF    (fn.CStr());

			}
			idx_out++;
			tout += dtOut;
			total_time = std::chrono::steady_clock::now() - start_whole;
			std::cout << "\nOutput No. " << idx_out << " at " << Time << " has been generated" << std::endl;
			std::cout << "Current Time Step = " <<deltat<<std::endl;
			std::cout << "Total time: "<<total_time.count() << ", Neigbour search time: " << clock_time_spent << ", Accel Calc time: " <<
			acc_time_spent <<
			std::endl;
			std::cout << "Max, Min, Avg temps: "<< Tmax << ", " << Tmin << ", " << (Tmax+Tmin)/2. <<std::endl;
			cout << "Max plastic strain: " <<max<< "in particle" << imax << endl;			
			double max_flux = 0.;
			for (size_t i=0; i<Particles.Size(); i++){
				if (Particles[i]->dTdt > max_flux)
					max_flux=Particles[i]->dTdt;
			}
			std::cout << "Max flux: "<< max_flux << std::endl;


      Vec3_t Max=1.e-10;
			int maxp[3];
			for (size_t i=0; i<Particles.Size(); i++) {
			if (Particles[i]->Displacement(0)*Particles[i]->Displacement(0) > Max(0)){Max(0) = Particles[i]->Displacement(0)*Particles[i]->Displacement(0); maxp[0]=i;}
				if (Particles[i]->Displacement(1)*Particles[i]->Displacement(0) > Max(1)){ Max(1) = Particles[i]->Displacement(1)*Particles[i]->Displacement(0);maxp[1]=i;}
				if (Particles[i]->Displacement(2)*Particles[i]->Displacement(0) > Max(2)){ Max(2) = Particles[i]->Displacement(2)*Particles[i]->Displacement(0);maxp[2]=i;}
			}
			cout << "Max Displacements: "<< sqrt(Max(0))<< ", "<<
																			sqrt(Max(1))<<", "<<
																			sqrt(Max(2))<<", "<<
																			endl;		
			//cout << "In particle" << maxp[0]<<", " << maxp[1]<<" , "<<maxp[2];
		}

		
	}//Main while
	

	std::cout << "\n--------------Solving is finished---------------------------------------------------" << std::endl;

}

}; //SPH