#include "Domain.h"

namespace SPH {

inline void Domain::CalcTempInc () {
	double di=0.0,dj=0.0,mi=0.0,mj=0.0;
	//#pragma omp parallel for schedule (static) num_threads(Nproc)
	for ( size_t k = 0; k < Nproc ; k++) {
		Particle *P1,*P2;
		Vec3_t xij;
		double h,GK;
		//TODO: DO THE LOCK PARALLEL THING
		// Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
		double temp=0;
		for (size_t a=0; a<SMPairs[k].Size();a++) {
			P1	= Particles[SMPairs[k][a].first];
			P2	= Particles[SMPairs[k][a].second];
			xij	= Particles[P1]->x-Particles[P2]->x;
			h	= (Particles[P1]->h+Particles[P2]->h)/2.0;
			GK	= GradKernel(Dimension, KernelType, norm(xij)/h, h);
			di = P1->Density; mi = P1->Mass;
			dj = P2->Density; mj = P2->Mass;
			//Frasier  Eqn 3.99 dTi/dt= 1/(rhoi_CPi) * Sum_j(mj/rho_j * 4*ki kj/ (ki + kj ) (Ti - Tj)  ) 
			temp += mj/dj * 4. * ( P1->k_T * P2->k_T) / (P1->k_T + P2->k_T) * ( P1->T * P2->T) * dot( xij , GK*xij );
			omp_set_lock(&P1->my_lock);
				P1->T		+= 1./(di*P1->cp_T) * ( temp + P1->q_conv );
			omp_unset_lock(&P1->my_lock);
			// Locking the particle 2 for updating the properties
			omp_set_lock(&P2->my_lock);
				P2->T		-= 1./(di*P2->cp_T) * ( temp + P2->q_conv );
			omp_unset_lock(&P2->my_lock);
		}
		// for (size_t i=0; i<FSMPairs[k].Size();i++)
			// CalcForce2233(Particles[FSMPairs[k][i].first],Particles[FSMPairs[k][i].second]);
	}//Nproc

}
inline void Domain::CalcConvHeat ( const int &id ){ //TODO: Detect Free Surface Elements
	double dS2;
	//Fraser Eq 3-121 
	#pragma omp parallel for schedule (static) num_threads(Nproc)
		for (size_t i=0; i<Particles.Size(); i++){	//Like in Domain::Move
			if ( Particles[i]->ID == id ) {
				dS2 = pow(Particles[i]->Mass/Particles[i]->Density,0.666666666);
				Particles[i]->q_conv=Particles[i]->Density * Particles[i]->h_conv * dS2 * (Particles[i]->T_inf - Particles[i]->T);
			}
		}	
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
	


	while (Time<tf && idx_out<=maxidx) {

		auto start_task = std::chrono::system_clock::now();
		clock_beg = clock();
		MainNeighbourSearch();
		clock_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		auto end_task = std::chrono::system_clock::now();
		 neighbour_time = /*std::chrono::duration_cast<std::chrono::seconds>*/ (end_task- start_task);
		//std::cout << "neighbour_time (chrono, clock): " << clock_time_spent << ", " << neighbour_time.count()<<std::endl;
		GeneralBefore(*this);
		clock_beg = clock();
		PrimaryComputeAcceleration();
		LastComputeAcceleration();
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
		}

		AdaptiveTimeStep();
		Move(deltat);
		Time += deltat;
		if (BC.InOutFlow>0) InFlowBCLeave(); else CheckParticleLeave ();
		CellReset();
		ListGenerate();
		
	}
	

	std::cout << "\n--------------Solving is finished---------------------------------------------------" << std::endl;

}

};
