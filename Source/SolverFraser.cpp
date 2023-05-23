namespace SPH {
  
///////////////////////////////////////
//This calculates all in same time step
//IMPORTANT: x is updated with previous accel!
// 1. x(t+dt) = x +v(t) dt + 1/2 a(t) dt^2
// 2. Calc a(t+dt) using x(t+dt)
// 3. v(t+dt) = v(t) dt + 1/2 (a(t) +a(t+dt)) * dt
///////////////////////////////////////
inline void Domain::SolveDiffUpdateFraser (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx) {
	std::cout << "\n--------------Solving---------------------------------------------------------------" << std::endl;

	size_t idx_out = 1;
	double tout = Time;

	//Initializing adaptive time step variables
	deltat = deltatint = deltatmin	= dt;

	clock_t clock_beg;
  double clock_time_spent,start_acc_time_spent, nb_time_spent ,pr_acc_time_spent,acc_time_spent, 
          contact_time_spent, trimesh_time_spent, bc_time_spent,
          mov_time_spent,stress_time_spent,energy_time_spent, dens_time_spent;
          
  clock_time_spent = contact_time_spent = acc_time_spent = stress_time_spent = energy_time_spent = dens_time_spent = mov_time_spent = 0.;	

	InitialChecks();
	CellInitiate();
	ListGenerate();
	PrintInput(TheFileKey);
	TimestepCheck();
	WholeVelocity();
	
  ofstream ofprop("Prop.csv", std::ios::out);
	//Initial model output
	if (TheFileKey!=NULL) {
		String fn;
		fn.Printf    ("%s_Initial", TheFileKey);
		WriteXDMF    (fn.CStr());
		std::cout << "\nInitial Condition has been generated\n" << std::endl;
	}
	
  //#ifdef NONLOCK_SUM
  InitReductionArraysOnce();
  //#endif
  //cout << "aref "<< Aref[0][0]<<endl;
  
	unsigned long steps=0;
	unsigned int first_step;
	
	int ts_i=0;

	bool isfirst = true;
	bool isyielding = false;
  
  //BEFORE CONTACT SURFACE SEARCH
	if (contact){
		for (int i=0; i<Particles.Size(); i++)
			Particles [i] -> ID_orig = Particles [i] -> ID;
	}

  double dS;
	if (contact) { //Calculate particle Stiffness
		for (int i=0; i<Particles.Size(); i++){
			double bulk = Particles[i]->Cs * Particles[i]->Cs *Particles[i]-> Density;  //RESTORE ORIGINAL BULK
			dS = pow(Particles[i]->Mass/Particles[i]->Density,0.33333); //Fraser 3-119
			Particles [i] -> cont_stiff = 9. * bulk * Particles [i]->G / (3. * bulk + Particles [i]->G) * dS;  //Fraser Thesis, Eqn. 3-153
		}		
    dS = pow(Particles[0]->Mass/Particles[0]->Density,0.33333);
		cout << "dS, psi_cont, Contact Stiffness" << dS << ", " 
    << Particles[0]->Cs * Particles [0] ->Mass / dS << ", " << Particles [0] -> cont_stiff <<endl;
		min_force_ts = deltat;
		MainNeighbourSearch();
    //CalcPairPosList();          //Only for h update
    //UpdateSmoothingLength();
		SaveNeighbourData();				//Necesary to calulate surface! Using Particle->Nb (count), could be included in search
		CalculateSurface(1);				//After Nb search	
	}
  
  ClearNbData();

        
	//Print history
	std::ofstream of("History.csv", std::ios::out);
  of << "Displacement, pl_strain, eff_strain_rate, sigma_eq, sigmay, contforcesum"<<endl;
  
  bool check_nb_every_time = false;

  cout << "Main Loop"<<endl;
  
  if (nonlock_sum && gradKernelCorr){
    cout << "WARNING: Nishimura summation is not working with Gradient Kernel Correction. Summation changed to Locking." <<endl;
    nonlock_sum = false;
  }
    
  int ct=30;
  std::chrono::duration<double> total_time;
	auto start_whole = std::chrono::steady_clock::now();  
  
  std::vector <Vec3_t> prev_acc(Particles.Size()); 
  #pragma omp parallel for schedule (static) num_threads(Nproc)	//LUCIANO//LIKE IN DOMAIN->MOVE
  for (int i=0; i<Particles.Size(); i++){
    prev_acc[i] = Particles[i]->a;
  }
  
  cout << std::setprecision(3)<< "Total allocated memory: " <<sizeof(Particle) * Particles.Size() * 1.0e-6 << " MB. "<<endl;
  
  while (Time<=tf && idx_out<=maxidx) {
  
		StartAcceleration(0.);

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

    max_disp = Vec3_t(0.,0.,0.);
		for (int i=0; i < solid_part_count; i++){
      for (int j=0;j<3;j++)
        if (Particles[i]->Displacement[j] * Particles[i]->Displacement[j]>max_disp[j]){
          max_disp[j] = Particles[i]->Displacement [j] * Particles[i]->Displacement [j];
          imax=i;
			}
		}
    for (int j=0;j<3;j++) max_disp[j] = sqrt(max_disp[j]);
    
    // // // ATTENTION! COULD BE LARGE DISPLACEMENTS AND SMALL STRAINS 
    // // // EXAMPLE COMPRESSION WITH NO FRICTION, SO CONTACTS NBs SHOULD BE RECALCULATED
    if (norm(max_disp) > 0.1 * hmax){
      if (!check_nb_every_time)
        cout << "Checking Nb Every step now."<<endl;
      check_nb_every_time = true;
    }
    else 
      check_nb_every_time = false;
		
		if (max > MIN_PS_FOR_NBSEARCH && !isyielding){ //First time yielding, data has not been cleared from first search
			ClearNbData(); 
			MainNeighbourSearch/*_Ext*/();
      CalcPairPosList();  
      SaveNeighbourData();
			if (contact) ContactNbUpdate(this);
			isyielding  = true ;
		}

		if ( max > MIN_PS_FOR_NBSEARCH || isfirst || check_nb_every_time){	//TO MODIFY: CHANGE
			if ( ts_i == 0 ){

				if (m_isNbDataCleared){
          clock_beg = clock();
					MainNeighbourSearch/*_Ext*/();
          //
          CalcPairPosList(); //For min TS Vel
          if (h_update){                
            UpdateSmoothingLength();          
          }
          //#ifdef NONLOCK_TEST 
          //CheckParticlePairs(0);
          //#endif
          SaveNeighbourData();
          //cout << "nb search"<<endl;
          nb_time_spent+=(double)(clock() - clock_beg) / CLOCKS_PER_SEC;
          
          if (gradKernelCorr){
            CalcGradCorrMatrix();	}
          
          if (contact) {
            clock_beg = clock();
            ContactNbUpdate(this);
            // if (isfirst)
              // CalcContactInitialGap(); //BEFORE! contactnb
            contact_time_spent +=(double)(clock() - clock_beg) / CLOCKS_PER_SEC;
          }
				}// ts_i == 0				
			}
		
    } //( max > MIN_PS_FOR_NBSEARCH || isfirst ){	//TO MODIFY: CHANGE

		// //NEW, gradient correction
			if (isfirst) {
				if (gradKernelCorr){
          CalcGradCorrMatrix();	}
				
			}		
		//std::cout << "neighbour_time (chrono, clock): " << clock_time_spent << ", " << neighbour_time.count()<<std::endl;

		GeneralBefore(*this);
		PrimaryComputeAcceleration();
    
    clock_beg = clock();
    //If density is calculated AFTER displacements, it fails
    CalcDensInc(); //TODO: USE SAME KERNEL?
    //#ifdef NONLOCK_SUM
    if (nonlock_sum) 
      DensReduction();
    //#endif
    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (size_t i=0; i<Particles.Size(); i++){
      //Particles[i]->UpdateDensity_Leapfrog(deltat);
      Particles[i]->Density += deltat*Particles[i]->dDensity;
    }    
    dens_time_spent+=(double)(clock() - clock_beg) / CLOCKS_PER_SEC;


		clock_beg = clock();
    CalcRateTensors();  //With v and xn+1

    //#ifdef NONLOCK_SUM
    if (nonlock_sum)
    RateTensorsReduction();
    //#endif
    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (size_t i=0; i<Particles.Size(); i++){
      //Particles[i]->Mat2Leapfrog(deltat); //Uses density  
      Particles[i]->CalcStressStrain(deltat); //Uses density  
    } 
    stress_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;		
    
    
    ////// 9. Update Density and Stress: already done
   
    //////////// 11 to 13.
    clock_beg = clock();
    //cout << "Particle 0 accel " << Particles[0]->a<<endl;
    CalcAccel(); //Nor density or neither strain rates
    //CalcAccelPP();
    //cout << "part 2000 acc "<<Particles[2000]->a<<endl;
    //#ifdef NONLOCK_SUM
    if (nonlock_sum) AccelReduction();
    else {
      if (dom_bid_type == AxiSymmetric)
      for (int i=0; i<solid_part_count;i++){
        Particles[i]->a *= 2.0 * M_PI;
        Particles[i]->a[0] -= Particles[i]->Sigma(2,2); //WANG Eqn. 40
      }    
    }
    //#endif
		acc_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
    
    cout << "calculate contact " <<endl;
    ///// 8. CONTACT FORCES
    clock_beg = clock(); 
    if (contact) CalcContactForcesWang();
    contact_time_spent +=(double)(clock() - clock_beg) / CLOCKS_PER_SEC;
    //if (contact) CalcContactForces2();
    cout << "done ."<<endl;
    
    //14. Add contact
    // if (contact ){
    // #pragma omp parallel for schedule (static) num_threads(Nproc)
    // for (size_t i=0; i<Particles.Size(); i++)
      // Particles[i]->a += Particles[i] -> contforce / Particles[i] -> Mass; 
    // }
    //// UPDATE VEL AND POS
    clock_beg = clock();  
    //BEFORE
    Vec3_t du;    
    GeneralAfter(*this);//Reinforce BC vel   
    //CorrectVelAcc();
    MoveGhost(); 
   #pragma omp parallel for schedule (static) private(du) num_threads(Nproc)
    for (size_t i=0; i<Particles.Size(); i++){
      Particles[i]->x_prev = Particles[i]->x;
      du = (Particles[i]->v + Particles[i]->VXSPH)*deltat + 0.5 * Particles[i]->a *deltat*deltat;
      Particles[i]->Displacement += du;
      Particles[i]->x += du;
    }
    mov_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;  
    
    
		
    double factor = 1.;

    clock_beg = clock();

    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (size_t i=0; i<Particles.Size(); i++){
      Particles[i]->v += Particles[i]->a * deltat;
    }
    //CorrectVelAcc();
    MoveGhost();   

    GeneralAfter(*this);//Reinforce BC vel   
    mov_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;  

    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (size_t i=0; i<Particles.Size(); i++)
      prev_acc[i] = Particles[i]->a;

		CalcPlasticWorkHeat(deltat);   //Before Thermal increment because it is used
    ThermalCalcs(deltat);

		clock_beg = clock();        
    CalcKinEnergyEqn();    
    CalcIntEnergyEqn();    
    energy_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
    
		steps++;
    if (ct == 30) ct = 0; else ct++;
    


		Time += deltat;
    clock_beg = clock();      
    if (contact){

      if (contact_mesh_auto_update) {
        for (int m=0; m<trimesh.size();m++)
          trimesh[m]->Update (deltat); //Update Node Pos, NOW includes PosCoeff and normals        
      }
      //cout << "Updating contact particles"<<endl;
      UpdateContactParticles(); //Updates normal and velocities
		}
    contact_time_spent +=(double)(clock() - clock_beg) / CLOCKS_PER_SEC;
		
		if (max>MIN_PS_FOR_NBSEARCH){	//TODO: CHANGE TO FIND NEIGHBOURS
			if ( ts_i == (ts_nb_inc - 1) ){
				ClearNbData();
			}

			ts_i ++;
			if ( ts_i > (ts_nb_inc - 1) ) 
				ts_i = 0;
		
		}
    
    if (Particles[0]->FirstStep)
    for (size_t i=0; i<Particles.Size(); i++){
      Particles[i]->FirstStep = false;
    }
		if (isfirst) isfirst = false;

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
			std::cout << "\n---------------------------------------\n Total CPU time: "<<total_time.count() << endl;
      double acc_time_spent_perc = acc_time_spent/total_time.count();
      std::cout << std::setprecision(2);
      cout << "Calculation Times\nAccel: "<<acc_time_spent_perc<<"%,  ";
      cout << "Density: "<<dens_time_spent/total_time.count()<<"%,  ";
      cout << "Stress: "  <<stress_time_spent/total_time.count()<<"%,  "<<endl;
      cout << "Energy: "  <<energy_time_spent/total_time.count()<<"%,  ";
      cout << "Contact: " <<contact_time_spent/total_time.count()<<"%,  ";
      cout << "Nb: "      <<nb_time_spent/total_time.count()<<"%,  ";
      cout << "Update: " <<mov_time_spent/total_time.count()<<"%,  ";
      cout <<endl;
      
			std::cout << "Output No. " << idx_out << " at " << Time << " has been generated" << std::endl;
			std::cout << "Current Time Step = " <<deltat<<std::endl;
			cout << "Max plastic strain: " <<max<< "in particle" << imax << endl;
      if (max > 0.)
        cout<<"Plastic Work "<<plastic_work<<endl;
      cout.precision(6);
			cout << "Max Displacements (No Cont Surf): "<<max_disp<<endl;
      if (contact) {
        cout<<"Contact Force Sum "<<contact_force_sum<<", Reaction Sum "<< contact_reaction_sum<<endl;
        cout<<"Contact Friction Work "<<contact_friction_work<<endl;
        cout<<"External Forces Work "<< ext_forces_work<<endl;
        if (cont_heat_cond)
          cout << "Total contact heat flux" << accum_cont_heat_cond <<endl;
      }

      cout << "Int Energy: " << int_energy_sum << ", Kin Energy: " << kin_energy_sum<<endl;
      if (thermal_solver)
        std::cout << "Max, Min, Avg temps: "<< m_maxT << ", " << m_minT << ", " << (m_maxT+m_minT)/2. <<std::endl;      
      
      ofprop <<getTime() << ", "<<m_scalar_prop<<endl;
			
      for (int p=0;p<Particles.Size();p++){
				if (Particles[p]->print_history)
          of << Particles[p]->Displacement << ", "<<Particles[p]->pl_strain<<", "<<Particles[p]->eff_strain_rate<<", "<< 
          Particles[p]->Sigma_eq<<", "  <<  Particles[p]->Sigmay << ", " <<
          contact_force_sum << endl;
			}
		}
    if (auto_ts)      CheckMinTSVel();
    if (auto_ts_acc)  CheckMinTSAccel();
    if (auto_ts || auto_ts_acc)  AdaptiveTimeStep();
    
	
	}


	of.close(); //History 
  ofprop.close(); //Scalar prop
	
	std::cout << "\n--------------Solving is finished---------------------------------------------------" << std::endl;

}

};