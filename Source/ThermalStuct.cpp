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
	double clock_time_spent,pr_acc_time_spent,acc_time_spent;
	double neigbour_time_spent_per_interval=0.;
	
	clock_time_spent=pr_acc_time_spent=acc_time_spent=0.;


	//Initial model output
	if (TheFileKey!=NULL) {
		String fn;
		fn.Printf    ("%s_Initial", TheFileKey);
		WriteXDMF    (fn.CStr());
		std::cout << "\nInitial Condition has been generated\n" << std::endl;
	}
	

	unsigned long steps=0;
	unsigned int first_step;
	
	int ts_nb_inc=5;	// Always > 0
	int ts_i=0;

	bool isfirst = true;
	bool isyielding = false;
	
	CalcConvHeat();
	CalcTempInc();
	
	while (Time<=tf && idx_out<=maxidx) {
		StartAcceleration(Gravity);
		//if (BC.InOutFlow>0) InFlowBCFresh();
		auto start_task = std::chrono::system_clock::now();
		

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

			MainNeighbourSearch/*_Ext*/();
			isyielding  = true ;
		}
		if ( max > MIN_PS_FOR_NBSEARCH || isfirst ){	//TO MODIFY: CHANGE
			if ( ts_i == 0 ){
				clock_beg = clock();
				if (m_isNbDataCleared)
					MainNeighbourSearch/*_Ext*/();

				neigbour_time_spent_per_interval += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
			}
			isfirst = false;
		}
		
		
		// for ( size_t k = 0; k < Nproc ; k++)		
			// cout << "Pares: " <<SMPairs[k].Size()<<endl;

		std::vector <int> nb(Particles.Size());
		std::vector <int> nbcount(Particles.Size());
		for ( size_t k = 0; k < Nproc ; k++) {
			for (size_t a=0; a<SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			//cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
				nb[SMPairs[k][a].first ]+=1;
				nb[SMPairs[k][a].second]+=1;
				
			}
		}	
			for (int p=0;p<Particles.Size();p++){
			Particles[p]->Nb=nb[p];
		}
	// for (int i=0;i<nb.size();i++)
		// cout << "Neigbour "<< i <<": "<<nb[i]<<endl;

		auto end_task = std::chrono::system_clock::now();
		 neighbour_time = /*std::chrono::duration_cast<std::chrono::seconds>*/ (end_task- start_task);
		//std::cout << "neighbour_time (chrono, clock): " << clock_time_spent << ", " << neighbour_time.count()<<std::endl;
		GeneralBefore(*this);
		clock_beg = clock();
		PrimaryComputeAcceleration();
		pr_acc_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		clock_beg = clock();
		LastComputeAcceleration();	//Define Strain Rate
		acc_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		GeneralAfter(*this);

		//Thermal calculations
		CalcConvHeat();
		CalcTempInc();	
		CalcThermalExpStrainRate();	//Add Thermal expansion Strain Rate Term
		
		for (size_t i=0; i<Particles.Size(); i++)
			Particles[i]->TempCalcLeapfrog(dt);
			
		steps++;
		//cout << "steps: "<<steps<<", time "<< Time<<", tout"<<tout<<endl;
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
			
			clock_time_spent += neigbour_time_spent_per_interval;
			std::cout << "Total CPU time: "<<total_time.count() << ", Neigbour search time: " << clock_time_spent << ", Pr Accel Calc time: " <<
			pr_acc_time_spent << "Las Acel Calc Time" << acc_time_spent<<
			std::endl;
						
			cout << "Max plastic strain: " <<max<< "in particle" << imax << endl;
			
			std::cout << "Steps count in this interval: "<<steps-first_step<<"Total Step count"<<steps<<endl;
			cout << "Total Nb search time in this interval: " << neigbour_time_spent_per_interval;
			cout << "Average Nb search time in this interval: " << neigbour_time_spent_per_interval/(float)(steps-first_step)<<endl;
			cout << "Avg Neighbour Count"<<AvgNeighbourCount()<<endl;
			first_step=steps;
			neigbour_time_spent_per_interval=0.;
		}
		
		if (auto_ts)
			AdaptiveTimeStep();
		Move(deltat);
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
	

	std::cout << "\n--------------Solving is finished---------------------------------------------------" << std::endl;

}


}; //SPH