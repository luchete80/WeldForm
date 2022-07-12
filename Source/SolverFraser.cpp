namespace SPH {
// VELOCITY AND POSITION UPDATE AT THE SAME TIME 
// THIS IS LIKE VERLET ALGORITHM 
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

}; //SPH