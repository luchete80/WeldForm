#include "Domain.h"
#include <vector>
#include "ThermalTest.cpp"

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
	
	#pragma omp parallel for schedule (static) num_threads(Nproc) //LUCIANO: THIS IS DONE SAME AS PrimaryComputeAcceleration
	for ( int k = 0; k < Nproc ; k++) {
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

      if (dom_bid_type == AxiSymmetric){ //CALCULATED DENSITY
        di/=(2.0*M_PI*P1->x(0));
        dj/=(2.0*M_PI*P2->x(0));
      }
			
			//Frasier  Eqn 3.99 dTi/dt= 1/(rhoi_CPi) * Sum_j(mj/rho_j * 4*ki kj/ (ki + kj ) (Ti - Tj)  ) 
			//LUCIANO: TODO EXCLUDE THIS PRODUCT
			double m, mc[2];
			if (gradKernelCorr){
				Mat3_t GKc[2];
				GKc[0] = P1->gradCorrM;
				GKc[1] = P2->gradCorrM;
				
				// if (SMPairs[k][a].first == 723){
					// cout << "Original GK * xij"<<GK * xij<<endl;
				// }
				//Left in vector form and multiply after??
				for (int i=0;i<2;i++){
					Vec3_t v;
					Mult (GKc[i], GK * xij, v);

					// if (i==0)
					// cout << P1->Nb<<endl;
					// else
					// cout << P2->Nb<<endl;
					mc[i]=mj/dj * 4. * ( P1->k_T * P2->k_T) / (P1->k_T + P2->k_T) * ( P1->T - P2->T) * dot( xij , v  )/ (norm(xij)*norm(xij));
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
	//Another test
	// for (int i=0; i<Particles.Size(); i++){
		// //double GK = temp[i] * GradKernel(Dimension, 0, 0., h);
		// Vec3_t GK_c; 
		// Mult(dom.Particles[i]->gradCorrM,temp[i],GK_c);
	// }
	
	//TODO: MULTIPLY CORRECTED GRADIENT HERE AFTER ALL SUM 
		//temp [i];
	
	double max = 0;
	int imax;
  
  double f;
  double frw, fr_temp = 0., cond_temp = 0.0;
  double plw, pl_sum = 0.;
  
	#pragma omp parallel for schedule (static) num_threads(Nproc)	private (f)//LUCIANO//LIKE IN DOMAIN->MOVE
  for (int i=0; i < solid_part_count; i++){
	//for (int i=0; i<Particles.Size(); i++){
		//cout << "temp "<<temp[i]<<endl;
    double d = Particles[i]->Density;
    if (dom_bid_type == AxiSymmetric) //CALCULATED DENSITY
      d/=(2.0*M_PI*Particles[i]->x(0));

      
		f = 1./(d * Particles[i]->cp_T ); //[ºC m^3/J]
    Particles[i]->dTdt = f * ( temp[i] + Particles[i]->q_conv + Particles[i]->q_source + Particles[i]->q_plheat * pl_work_heat_frac + Particles[i]->q_cont_conv);	
    
    plw = f * Particles[i]->q_plheat;
    pl_sum += plw;
    // Particles[i]->dTdt += plw;
    
		if (contact){
      if (cont_heat_fric) {
        frw = f * Particles[i]->q_fric_work; //[ºC m^3/J] x J/[s m3] = ºC/s
    
        //omp_set_lock(&dom_lock); 
        // #pragma omp atomic
        // fr_temp += Particles[i]->q_fric_work * Particles[i]->Mass / Particles[i]->Density; //TODO: CHECK  -- J/[s m3] x m3
        //omp_unset_lock(&dom_lock);	
        
        omp_set_lock(&Particles[i]->my_lock);
        Particles[i]->dTdt += frw; //[J/(kg.s)] / [J/(kg.K)]]
        omp_unset_lock(&Particles[i]->my_lock);
      }
    }
		if (Particles[i]->dTdt > max){
			max= Particles[i]->dTdt;
			imax=i;
		}
	}
  
  
  if (contact){

    //REAL VOLUME IN AXISYMM is 2PI*r
    for (int i=0; i < solid_part_count; i++){
      double d = Particles[i]->Density;
      if (dom_bid_type == AxiSymmetric) //CALCULATED DENSITY
        d/=(2.0*M_PI*Particles[i]->x(0));
        
      fr_temp += Particles[i]->q_fric_work * Particles[i]->Mass / d;
    if (cont_heat_cond)
      cond_temp += Particles[i]->q_cont_conv * Particles[i]->Mass / d;
    }
  }
  for (int i=solid_part_count; i<Particles.Size(); i++){
    Particles[i]->dTdt = 0.;
  }
  contact_friction_work += fr_temp * deltat;
  accum_cont_heat_cond  += cond_temp * deltat;
  
  //plastic_work += pl_sum * deltat;
  
	//cout << "Max dTdt: " << max <<"in particle: " << imax<<endl;
	
}

inline void Domain::CalcTempIncSOA () {
	double di=0.0,dj=0.0,mi=0.0,mj=0.0;
	
	std::vector < double> temp(Particles.Size());
	
	//cout << "calc temp inc"<<endl;
	#pragma omp parallel for schedule (static) num_threads(Nproc) //LUCIANO: THIS IS DONE SAME AS PrimaryComputeAcceleration
	for ( int k = 0; k < Nproc ; k++) {
		int P1,P2;
		Vec3_t xij;
		double h,GK;
		//TODO: DO THE LOCK PARALLEL THING
		// Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
		//std::vector <double> temp()=0;
		//cout << "fixed pair size: "<<FSMPairs[k].Size()<<endl;
		//cout << "Particles size: " << Particles.Size()<<endl;
		for (size_t a=0; a<SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			//cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
			P1	= SMPairs[k][a].first;
			P2	= SMPairs[k][a].second;
			xij	= *m_x[P1] - (*m_x[P2]);
			h	= (*m_h[P1] + (*m_h[P2]))/2.0;
			GK	= GradKernel(Dimension, KernelType, norm(xij)/h, h);	
			
			
			di = *m_rho[P1]; mi = *m_mass[P1];
			dj = *m_rho[P2]; mj = *m_mass[P2];
			
			//cout << "calc"<<endl;
			//Frasier  Eqn 3.99 dTi/dt= 1/(rhoi_CPi) * Sum_j(mj/rho_j * 4*ki kj/ (ki + kj ) (Ti - Tj)  ) 
			//LUCIANO: TODO EXCLUDE THIS PRODUCT
			double m, mc[2];
			if (gradKernelCorr){
				// Mat3_t GKc[2];
				// GKc[0] = P1->gradCorrM;
				// GKc[1] = P2->gradCorrM;
				
				// // if (SMPairs[k][a].first == 723){
					// // cout << "Original GK * xij"<<GK * xij<<endl;
				// // }
				// //Left in vector form and multiply after??
				// for (int i=0;i<2;i++){
					// Vec3_t v;
					// Mult (GKc[i], GK * xij, v);
					// // if (SMPairs[k][a].first == 723)
					// // cout << "Orig, Corr GK * xij, Nb"<<GK * xij<<", "<< v;
					// // if (i==0)
					// // cout << P1->Nb<<endl;
					// // else
					// // cout << P2->Nb<<endl;
					// mc[i]=mj/dj * 4. * ( P1->k_T * P2->k_T) / (P1->k_T + P2->k_T) * ( P1->T - P2->T) * dot( xij , v  )/ (norm(xij)*norm(xij));
				// }				
			} else {
				m = mj/dj * 4. * ( *m_kT[P1] * (*m_kT[P2])) / (*m_kT[P1] + *m_kT[P2]) * ( *m_T[P1] - (*m_T[P2])) * dot( xij , GK*xij )/ (norm(xij)*norm(xij));
				mc[0]=mc[1]=m;
			}
			//omp_set_lock(&P1->my_lock);
			temp [P1]  += mc[0];
			temp [P2] -= mc[1];
		}
	}//Nproc
	//Another test
	// for (int i=0; i<Particles.Size(); i++){
		// //double GK = temp[i] * GradKernel(Dimension, 0, 0., h);
		// Vec3_t GK_c; 
		// Mult(dom.Particles[i]->gradCorrM,temp[i],GK_c);
	// }
	
	//TODO: MULTIPLY CORRECTED GRADIENT HERE AFTER ALL SUM 
		//temp [i];
	
	//cout << "end temp calculating dTdt"<<endl;
	double max = 0;
	int imax;
	#pragma omp parallel for schedule (static) num_threads(Nproc)	//LUCIANO//LIKE IN DOMAIN->MOVE
	for (int i=0; i<Particles.Size(); i++){
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
			Particles[i]->q_conv = Particles[i]->Density * Particles[i]->h_conv * dS2 * (Particles[i]->T_inf - Particles[i]->T)/Particles[i]->Mass;
			
			if (Particles[i]->q_conv>max){
				max= Particles[i]->q_conv;
				imax=i;
			}
			//cout << "Particle  conv"<<Particles[i]->q_conv<<endl;
		}
	}		
	//cout << "Max Convection: " << max <<"in particle " << imax <<endl;
	//cout << "Applied convection to "<< i << " Particles"<<endl;
}

inline void Domain::CalcConvHeatSOA (){ //TODO: Detect Free Surface Elements
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
			dS2 = pow(*m_mass[i]/(*m_rho[i]),0.666666666);
			//cout << "dS2" <<dS2<<endl;
			//cout << Particles[i]->Density<<endl;
			//Fraser Eq 3.121
			*m_qconvT[i] = *m_rho[i] * (*m_hcT[i]) * dS2 * (*m_Tinf[i] - *m_T[i])/(*m_mass[i]);
			
			if (*m_qconvT[i] > max){
				max 	= *m_qconvT[i];
				imax	=	i;
			}
			//cout << "Particle  "<<Particles[i]->Mass<<endl;
		}
	}		
	//cout << "Max Convection: " << max <<"in particle " << imax <<endl;
	//cout << "Applied convection to "<< i << " Particles"<<endl;
}

//PREVIOUS THERMAL SOLVER
// inline void Domain::ThermalSolve (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx) {
	// std::cout << "\n--------------Solving---------------------------------------------------------------" << std::endl;

	// size_t idx_out = 1;
	// double tout = Time;

	// //Initializing adaptive time step variables
	// deltat = deltatint = deltatmin	= dt;
	
	// auto start_whole = std::chrono::steady_clock::now();

	// InitialChecks();
	// CellInitiate();
	// ListGenerate();
	// PrintInput(TheFileKey);
	// TimestepCheck();
	// WholeVelocity();
	
	// std::chrono::duration<double> total_time,neighbour_time;
	
	// clock_t clock_beg;
	// double clock_time_spent,acc_time_spent;
	
	// clock_time_spent=acc_time_spent=0.;


	// //Initial model output
	// if (TheFileKey!=NULL) {
		// String fn;
		// fn.Printf    ("%s_Initial", TheFileKey);
		// WriteXDMF    (fn.CStr());
		// std::cout << "\nInitial Condition has been generated\n" << std::endl;
	// }
	
	// MainNeighbourSearch();
	// SaveNeighbourData();
	// cout << "Avg Nb Count: "<<AvgNeighbourCount()<<endl;
	
				// String fn;
				// fn.Printf    ("%s_%04d", TheFileKey, idx_out);
				// WriteXDMF    (fn.CStr());
	// if (gradKernelCorr)
		// CalcGradCorrMatrix();
	// for ( size_t k = 0; k < Nproc ; k++) 
		// cout << "Pares: " <<SMPairs[k].Size()<<endl;

	// // ONLY FOR TEST
  // cout << "Init red arrays"<<endl;
  // InitReductionArraysOnce();
  // CalcPairPosList();  
  // cout << "Done."<<endl;
  
	// cout << "Calc conv "<<endl;
	// //CalcConvHeatSOA();
	// cout << "Done. "<<endl;
	// //CalcTempIncSOA();
  // //CalcTempInc();
  // CalcTempIncPP();
  // //CalcConvHeat();
	// cout << "End."<<endl;	
	// while (Time<tf && idx_out<=maxidx) {

		// auto start_task = std::chrono::system_clock::now();
		// clock_beg = clock();
		// clock_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		// auto end_task = std::chrono::system_clock::now();
		 // neighbour_time = /*std::chrono::duration_cast<std::chrono::seconds>*/ (end_task- start_task);
		// //std::cout << "neighbour_time (chrono, clock): " << clock_time_spent << ", " << neighbour_time.count()<<std::endl;
		// clock_beg = clock();

		// // CalcConvHeat();
		// // CalcTempInc();
		// //TODO Add 
		// double max=0,min=1000.; 
    
    // #pragma omp parallel for schedule (static) num_threads(Nproc)
		// for (int i=0; i<Particles.Size(); i++){
			// Particles[i]->T+= dt*Particles[i]->dTdt;
			// //Particles[i]->TempCalcLeapfrog(dt);
			// //*m_T[i]+= (*m_dTdt[i])*dt;
			// if (Particles[i]->T > max)
				// max=Particles[i]->T;
			// if (Particles[i]->T < min)
				// min=Particles[i]->T;

			// // if (*m_T[i] > max)
				// // max = *m_T[i];
			// // if (*m_T[i] < min)
				// // min = *m_T[i];

		// }
		// // std::cout << "Max temp: "<< max << std::endl;

			
		// acc_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		// GeneralAfter(*this);
		// // output
		// if (Time>=tout){
			// if (TheFileKey!=NULL) {
				// String fn;
				// fn.Printf    ("%s_%04d", TheFileKey, idx_out);
				// WriteXDMF    (fn.CStr());

			// }
			// idx_out++;
			// tout += dtOut;
			// total_time = std::chrono::steady_clock::now() - start_whole;
			// std::cout << "\nOutput No. " << idx_out << " at " << Time << " has been generated" << std::endl;
			// std::cout << "Current Time Step = " <<deltat<<std::endl;
			// std::cout << "Total time: "<<total_time.count() << ", Neigbour search time: " << clock_time_spent << ", Accel Calc time: " <<
			// acc_time_spent <<
			// std::endl;
			// std::cout << "Max, Min, Avg temps: "<< max << ", " << min << ", " << (max+min)/2. <<std::endl;

    // cout << "Particle 0 dTdt "<<Particles[0]->dTdt<<endl;
    
			// double max_flux = 0.;
			// for (size_t i=0; i<Particles.Size(); i++){
				// if (Particles[i]->dTdt > max_flux)
					// max_flux=Particles[i]->dTdt;
			// }
			// std::cout << "Max flux: "<< max_flux << std::endl;

		// }

		// //AdaptiveTimeStep();
		
		// //CalcConvHeatSOA();
		// CalcTempInc();
    // //CalcConvHeat();
    // //CalcTempIncSOA();

		// Time += deltat;
		
	// }
	

	// std::cout << "\n--------------Solving is finished---------------------------------------------------" << std::endl;

// }



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
	SaveNeighbourData();
	cout << "Avg Nb Count: "<<AvgNeighbourCount()<<endl;
	
				String fn;
				fn.Printf    ("%s_%04d", TheFileKey, idx_out);
				WriteXDMF    (fn.CStr());
	if (gradKernelCorr)
		CalcGradCorrMatrix();
	for ( size_t k = 0; k < Nproc ; k++) 
		cout << "Pares: " <<SMPairs[k].Size()<<endl;

	// ONLY FOR TEST
  cout << "Init red arrays"<<endl;
  InitReductionArraysOnce();
  CalcPairPosList();  
  cout << "Done."<<endl;
  
	cout << "Calc conv "<<endl;
	//CalcConvHeatSOA();
	cout << "Done. "<<endl;
	//CalcTempIncSOA();
  //CalcTempInc();
  CalcTempIncPP();
  //CalcConvHeat();
	cout << "End."<<endl;	
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
    
    // #pragma omp parallel for schedule (static) num_threads(Nproc)
		// for (int i=0; i<Particles.Size(); i++){
			// Particles[i]->T+= dt*Particles[i]->dTdt;
			// //Particles[i]->TempCalcLeapfrog(dt);
			// //*m_T[i]+= (*m_dTdt[i])*dt;
			// if (Particles[i]->T > max)
				// max=Particles[i]->T;
			// if (Particles[i]->T < min)
				// min=Particles[i]->T;

			// // if (*m_T[i] > max)
				// // max = *m_T[i];
			// // if (*m_T[i] < min)
				// // min = *m_T[i];

		// }
		// // std::cout << "Max temp: "<< max << std::endl;
    ThermalCalcs(dt);

    #pragma omp parallel for schedule (static) num_threads(Nproc)
		for (int i=0; i<Particles.Size(); i++){

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

    cout << "Particle 0 dTdt "<<Particles[0]->dTdt<<endl;
    
			double max_flux = 0.;
			for (size_t i=0; i<Particles.Size(); i++){
				if (Particles[i]->dTdt * Particles[i]->dTdt > max_flux)
					max_flux=Particles[i]->dTdt*Particles[i]->dTdt;
			}
			std::cout << "Max flux: "<< sqrt(max_flux) << std::endl;

		}

		//AdaptiveTimeStep();
		
		//CalcConvHeatSOA();
		//CalcTempInc();
    //CalcConvHeat();
    //CalcTempIncSOA();

		Time += deltat;
		
	}
	

	std::cout << "\n--------------Solving is finished---------------------------------------------------" << std::endl;

}


inline void Domain::ThermalCalcs(const double &dt){
  if (thermal_solver){
    m_maxT = 0.;
    m_minT =1000.; 

    #pragma omp parallel for schedule (static) num_threads(Nproc)    
    for (int i=0; i < solid_part_count; i++){
			Particles[i]->T += dt*Particles[i]->dTdt;
			//Particles[i]->TempCalcLeapfrog(dt);
			if (Particles[i]->T > m_maxT)
				m_maxT=Particles[i]->T;
			if (Particles[i]->T < m_minT)
				m_minT=Particles[i]->T;      
    }
    
			CalcConvHeat();
			CalcTempInc();
			CalcThermalExpStrainRate();	//Add Thermal expansion Strain Rate Term	
  
    if (contact && cont_heat_cond){
      #pragma omp parallel for schedule (static) num_threads(Nproc)    
      for (int i = solid_part_count; i < Particles.Size(); i++){
        Particles[i]->T += tot_cont_heat_cond[Particles[i]->mesh] / Particles[i]->mcp_t * deltat;
      }
    }
  
  }
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
