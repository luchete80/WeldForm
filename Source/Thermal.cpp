#include "Domain.h"
#include <vector>
using namespace std;

namespace SPH {

// inline void Domain::CalcTempInc () {
	// double di=0.0,dj=0.0,mi=0.0,mj=0.0;
	
	// std::vector < double> temp(Particles.Size());
	
	// #pragma omp parallel for schedule (static) num_threads(Nproc) //LUCIANO: THIS IS DONE SAME AS PrimaryComputeAcceleration
	// for ( size_t k = 0; k < Nproc ; k++) {
		// Particle *P1,*P2;
		// Vec3_t xij;
		// double h,GK;
		// //TODO: DO THE LOCK PARALLEL THING
		// // Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
		// //std::vector <double> temp()=0;
		// //cout << "fixed pair size: "<<FSMPairs[k].Size()<<endl;
		// //cout << "Particles size: " << Particles.Size()<<endl;
		// for (size_t a=0; a<SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			// //cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
			// P1	= Particles[SMPairs[k][a].first];
			// P2	= Particles[SMPairs[k][a].second];
			// xij	= P1->x - P2->x;
			// h	= (P1->h+P2->h)/2.0;
			// GK	= GradKernel(Dimension, KernelType, norm(xij)/h, h);
			// di = P1->Density; mi = P1->Mass;
			// dj = P2->Density; mj = P2->Mass;

			// //Frasier  Eqn 3.99 dTi/dt= 1/(rhoi_CPi) * Sum_j(mj/rho_j * 4*ki kj/ (ki + kj ) (Ti - Tj)  ) 
			// //LUCIANO: TODO EXCLUDE THIS PRODUCT
			// temp [SMPairs[k][a].first] += mj/dj * 4. * ( P1->k_T * P2->k_T) / (P1->k_T + P2->k_T) * ( P1->T - P2->T) * dot( xij , GK*xij );
			// omp_set_lock(&P1->my_lock);
				// P1->dTdt		+= 1./(di*P1->cp_T) * ( temp + P1->q_conv );
			// omp_unset_lock(&P1->my_lock);
			// // Locking the particle 2 for updating the properties
			// omp_set_lock(&P2->my_lock);
				// P2->dTdt		-= 1./(di*P2->cp_T) * ( temp + P2->q_conv );
			// omp_unset_lock(&P2->my_lock);
				// //cout << "temp: "<<dot( xij , GK*xij )<<endl;
		// }

	// }//Nproc

// }


//NON PARALLEL TEST VERSION
inline void Domain::CalcTempInc () {
	double di=0.0,dj=0.0,mi=0.0,mj=0.0;
	
	std::vector < double> temp(Particles.Size());
	
	//#pragma omp parallel for schedule (static) num_threads(Nproc) //LUCIANO: THIS IS DONE SAME AS PrimaryComputeAcceleration
	for ( size_t k = 0; k < Nproc ; k++) {
		Particle *P1,*P2;
		Vec3_t xij;
		double h,GK;
		//TODO: DO THE LOCK PARALLEL THING
		// Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
		//std::vector <double> temp()=0;
		//cout << "fixed pair size: "<<FSMPairs[k].Size()<<endl;
		//cout << "Particles size: " << Particles.Size()<<endl;
		for (size_t a=0; a<SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			//cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
			P1	= Particles[SMPairs[k][a].first];
			P2	= Particles[SMPairs[k][a].second];
			xij	= P1->x - P2->x;
			h	= (P1->h+P2->h)/2.0;
			GK	= GradKernel(Dimension, KernelType, norm(xij)/h, h);	
			
			
			di = P1->Density; mi = P1->Mass;
			dj = P2->Density; mj = P2->Mass;
			
			//Frasier  Eqn 3.99 dTi/dt= 1/(rhoi_CPi) * Sum_j(mj/rho_j * 4*ki kj/ (ki + kj ) (Ti - Tj)  ) 
			//LUCIANO: TODO EXCLUDE THIS PRODUCT
			double m, mc[2];
			if (gradKernelCorr){
				Mat3_t GKc[2];
				GKc[0] = GK * P1->gradCorrM;
				GKc[1] = GK * P2->gradCorrM;
				
				//Left in vector form and multiply after??
				for (int i=0;i<2;i++){
					Vec3_t v;
					Mult (GKc[i], xij, v);
					mc[i]=mj/dj * 4. * ( P1->k_T * P2->k_T) / (P1->k_T + P2->k_T) * ( P1->T - P2->T) * dot( xij , v )/ (norm(xij)*norm(xij));
				}				
			} else {
				m = mj/dj * 4. * ( P1->k_T * P2->k_T) / (P1->k_T + P2->k_T) * ( P1->T - P2->T) * dot( xij , GK*xij )/ (norm(xij)*norm(xij));
				mc[0]=mc[1]=m;
			}
			//omp_set_lock(&P1->my_lock);
			temp [SMPairs[k][a].first]  += mc[0];
			temp [SMPairs[k][a].second] -= mc[1];
		}
	}//Nproc
	
	//TODO: MULTIPLY CORRECTED GRADIENT HERE AFTER ALL SUM 
		//temp [i];
	
	double max = 0;
	int imax;
	#pragma omp parallel for schedule (static) num_threads(Nproc)	//LUCIANO//LIKE IN DOMAIN->MOVE
	for (int i=0; i<Particles.Size(); i++){
		//cout << "temp "<<temp[i]<<endl;
		Particles[i]->dTdt = 1./(Particles[i]->Density * Particles[i]->cp_T ) * ( temp[i] + Particles[i]->q_conv + Particles[i]->q_source);	
		if (Particles[i]->dTdt>max){
			max= Particles[i]->dTdt;
			imax=i;
		}
	}
	//cout << "Max dTdt: " << max <<"in particle: " << imax<<endl;
	
}

inline void Domain::CalcConvHeat (){ //TODO: Detect Free Surface Elements
	double dS2;
	//Fraser Eq 3-121 
	double max=0.;
	int imax;
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	
	#ifdef __GNUC__
	for (size_t i=0; i<Particles.Size(); i++){	//Like in Domain::Move
	#else
	for (int i=0; i<Particles.Size(); i++){//Like in Domain::Move
	#endif

	
		if ( Particles[i]->Thermal_BC==TH_BC_CONVECTION) {
			dS2 = pow(Particles[i]->Mass/Particles[i]->Density,0.666666666);
			//cout << "dS2" <<dS2<<endl;
			//cout << Particles[i]->Density<<endl;
			//Fraser Eq 3.121
			Particles[i]->q_conv=Particles[i]->Density * Particles[i]->h_conv * dS2 * (Particles[i]->T_inf - Particles[i]->T)/Particles[i]->Mass;
			if (Particles[i]->q_conv>max){
				max= Particles[i]->q_conv;
				imax=i;
			}
			//cout << "Particle  "<<Particles[i]->Mass<<endl;
		}
	}		
	//cout << "Max Convection: " << max <<"in particle " << imax <<endl;
	//cout << "Applied convection to "<< i << " Particles"<<endl;
}

inline void Domain::CalcPlasticWorkHeat (){ //TODO: Detect Free Surface Elements

	//Fraser Eq 3-106
	double max=0.;
	int imax;
	#pragma omp parallel for schedule (static) num_threads(Nproc)

	#ifdef __GNUC__
	for (size_t i=0; i<Particles.Size(); i++){	//Like in Domain::Move
	#else
	for (int i=0; i<Particles.Size(); i++){//Like in Domain::Move
	#endif
	
			// //cout << "dS2" <<dS2<<endl;
			// //cout << Particles[i]->Density<<endl;
			// Particles[i]->q_plheat=
			
					// double J2	= 0.5*(ShearStress(0,0)*ShearStress(0,0) + 2.0*ShearStress(0,1)*ShearStress(1,0) +
						// 2.0*ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,1)*ShearStress(1,1) +
						// 2.0*ShearStress(1,2)*ShearStress(2,1) + ShearStress(2,2)*ShearStress(2,2));
						
						// ;
			// if (Particles[i]->q_conv>max){
				// max= Particles[i]->q_conv;
				// imax=i;
			// }
			// //cout << "Particle  "<<Particles[i]->Mass<<endl;
		
	}		
	//cout << "Max Convection: " << max <<"in particle " << imax <<endl;
	//cout << "Applied convection to "<< i << " Particles"<<endl;
}

inline void Domain::ThermalSolve (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx) {
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
	double clock_time_spent,acc_time_spent;
	
	clock_time_spent=acc_time_spent=0.;


	//Initial model output
	if (TheFileKey!=NULL) {
		String fn;
		fn.Printf    ("%s_Initial", TheFileKey);
		WriteXDMF    (fn.CStr());
		std::cout << "\nInitial Condition has been generated\n" << std::endl;
	}
	
	MainNeighbourSearch();
	if (gradKernelCorr)
		CalcGradCorrMatrix();
	for ( size_t k = 0; k < Nproc ; k++) 
		cout << "Pares: " <<SMPairs[k].Size()<<endl;

	std::vector <int> nb(Particles.Size());
	std::vector <int> nbcount(Particles.Size());
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for ( int k = 0; k < Nproc ; k++) {
		for (int a=0; a<SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			nb[SMPairs[k][a].first ]+=1;
			nb[SMPairs[k][a].second]+=1;
		}
	}	
	for (int p=0;p<Particles.Size();p++){
		Particles[p]->Nb=nb[p];
	}

	unsigned long avg=0;
	for (int i=0;i<nb.size();i++) {
		avg+=nb[i];
	}
	avg/=nb.size();
	cout << "Avg Neigbour : "<<avg<<endl;
	
	CalcConvHeat();
	CalcTempInc();
	
	while (Time<tf && idx_out<=maxidx) {

		auto start_task = std::chrono::system_clock::now();
		clock_beg = clock();
		clock_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		auto end_task = std::chrono::system_clock::now();
		 neighbour_time = /*std::chrono::duration_cast<std::chrono::seconds>*/ (end_task- start_task);
		//std::cout << "neighbour_time (chrono, clock): " << clock_time_spent << ", " << neighbour_time.count()<<std::endl;
		clock_beg = clock();

		// CalcConvHeat();
		// CalcTempInc();
		//TODO Add 
		double max=0,min=1000.;
		for (size_t i=0; i<Particles.Size(); i++){
			//Particles[i]->T+= dt*Particles[i]->dTdt;
			Particles[i]->TempCalcLeapfrog(dt);
			if (Particles[i]->T > max)
				max=Particles[i]->T;
			if (Particles[i]->T < min)
				min=Particles[i]->T;
		}
		// std::cout << "Max temp: "<< max << std::endl;

			
		acc_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		GeneralAfter(*this);
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
			std::cout << "Max, Min, Avg temps: "<< max << ", " << min << ", " << (max+min)/2. <<std::endl;
			
			double max_flux = 0.;
			for (size_t i=0; i<Particles.Size(); i++){
				if (Particles[i]->dTdt > max_flux)
					max_flux=Particles[i]->dTdt;
			}
			std::cout << "Max flux: "<< max_flux << std::endl;

		}

		//AdaptiveTimeStep();
		
				CalcConvHeat();
		CalcTempInc();

		Time += deltat;
		
	}
	

	std::cout << "\n--------------Solving is finished---------------------------------------------------" << std::endl;

}

inline void Domain::ThermalSolve_wo_init (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx) {
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
	double clock_time_spent,acc_time_spent;
	
	clock_time_spent=acc_time_spent=0.;


	//Initial model output
	if (TheFileKey!=NULL) {
		String fn;
		fn.Printf    ("%s_Initial", TheFileKey);
		WriteXDMF    (fn.CStr());
		std::cout << "\nInitial Condition has been generated\n" << std::endl;
	}
	
	//MainNeighbourSearch();
	std::vector <int> nb(Particles.Size());
	std::vector <int> nbcount(Particles.Size());
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	for ( int k = 0; k < Nproc ; k++) {
		for (int a=0; a<SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			nb[SMPairs[k][a].first ]+=1;
			nb[SMPairs[k][a].second]+=1;
		}
	}	
	for (int p=0;p<Particles.Size();p++){
		Particles[p]->Nb=nb[p];
	}
	
		
	while (Time<tf && idx_out<=maxidx) {

		auto start_task = std::chrono::system_clock::now();
		clock_beg = clock();
		clock_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		auto end_task = std::chrono::system_clock::now();
		 neighbour_time = /*std::chrono::duration_cast<std::chrono::seconds>*/ (end_task- start_task);
		//std::cout << "neighbour_time (chrono, clock): " << clock_time_spent << ", " << neighbour_time.count()<<std::endl;
		clock_beg = clock();

		CalcConvHeat();
		CalcTempInc();
		//TODO Add 
		double max=0,min=1000.;
		for (size_t i=0; i<Particles.Size(); i++){
			Particles[i]->T+= dt*Particles[i]->dTdt;
			if (Particles[i]->T > max)
				max=Particles[i]->T;
			if (Particles[i]->T < min)
				min=Particles[i]->T;
		}
		// std::cout << "Max temp: "<< max << std::endl;

			
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
			std::cout << "Max, Min, Avg temps: "<< max << ", " << min << ", " << (max+min)/2. <<std::endl;
			
			double max_flux = 0.;
			for (size_t i=0; i<Particles.Size(); i++){
				if (Particles[i]->dTdt > max_flux)
					max_flux=Particles[i]->dTdt;
			}
			std::cout << "Max flux: "<< max_flux << std::endl;

		}

		AdaptiveTimeStep();

		Time += deltat;
		
	}
	

	std::cout << "\n--------------Solving is finished---------------------------------------------------" << std::endl;

}

};
