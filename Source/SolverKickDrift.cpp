////////////////////////////////////////////////////////
////////////////////// BASED ON RANDLES AND LIBERSKY(1996):
/////////////////////////////////////////////////////////
/// Randles and Libersky calculate density from current velocity, here is from the velocity at t+1/2
    // // // 1 CalcAccel(); //Nor density or neither strain rates
    // // // if (nonlock_sum)AccelReduction();

    // // // if (contact) CalcContactForcesWang();

    // // // 3. Particles[i]->v += Particles[i]->a*deltat/2.*factor;

    // // // 4. //If density is calculated AFTER displacements, it fails
    // // // CalcDensInc(); //TODO: USE SAME KERNEL?
    // // // Particles[i]->Density += deltat*Particles[i]->dDensity*factor;
    // // // 5. x += (Particles[i]->v + Particles[i]->VXSPH)*deltat*factor;
    
    // // // 6. Particles[i]->v += Particles[i]->a*deltat/2.*factor;
    // // // 7. CalcRateTensors();  //With v and xn+1
    // // // 8. Particles[i]->CalcStressStrain(deltat); //Uses density  


namespace SPH {
inline void Domain::SolveDiffUpdateKickDrift (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx) {
	std::cout << "\n--------------Solving---------------------------------------------------------------" << std::endl;

	size_t idx_out = 1;
	double tout = Time;

	//Initializing adaptive time step variables
	deltat = deltatint = deltatmin	= dt;

	clock_t clock_beg;
  double clock_time_spent,start_acc_time_spent, nb_time_spent ,pr_acc_time_spent,acc_time_spent, 
          contact_time_spent, trimesh_time_spent, bc_time_spent,
          mov_time_spent,stress_time_spent,energy_time_spent, dens_time_spent, thermal_time_spent,
          ini_time_spent;
  
          
  clock_time_spent = contact_time_spent = acc_time_spent = stress_time_spent = energy_time_spent = dens_time_spent = mov_time_spent = 0.;	
  double last_output_time;
  
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
  

  int ct=30;
  std::chrono::duration<double> total_time;
	auto start_whole = std::chrono::steady_clock::now();  
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

    Vec3_t max_disp = Vec3_t(0.,0.,0.);
		for (int i=0; i<Particles.Size(); i++){
      for (int j=0;j<3;j++)
        if (Particles[i]->Displacement[j]>max_disp[j]){
          max_disp[j] = Particles[i]->Displacement [j];
          imax=i;
			}
		}
    
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
          if (contact) {
            clock_beg = clock();
            ContactNbUpdate(this);
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
    //cout << "Particle 0 accel " << Particles[0]->a<<endl;
    CalcAccel(); //Nor density or neither strain rates
    //CalcAccelPP();
    //cout << "part 2000 acc "<<Particles[2000]->a<<endl;
    //#ifdef NONLOCK_SUM
    if (nonlock_sum)AccelReduction();
    //#endif
		acc_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
    GeneralAfter(*this); //Fix free accel
    
    clock_beg = clock(); 
    if (contact) CalcContactForcesWang();
    contact_time_spent +=(double)(clock() - clock_beg) / CLOCKS_PER_SEC;
    //if (contact) CalcContactForces2();
		
    double factor = 1.;
    // if (ct==30) factor = 1.;
    // else        factor = 2.;
    clock_beg = clock();
    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (int i=0; i<Particles.Size(); i++){
      Particles[i]->v += Particles[i]->a*deltat/2.*factor;
      //Particles[i]->LimitVel();
    }
    MoveGhost();   
    GeneralAfter(*this);//Reinforce BC vel   
    mov_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;  
    
    clock_beg = clock();
    //If density is calculated AFTER displacements, it fails
    CalcDensInc(); //TODO: USE SAME KERNEL?
    //#ifdef NONLOCK_SUM
    if (nonlock_sum) 
      DensReduction();
    //#endif    
    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (int i=0; i<Particles.Size(); i++){
      //Particles[i]->UpdateDensity_Leapfrog(deltat);
      Particles[i]->Density += deltat*Particles[i]->dDensity*factor;
    }    
    dens_time_spent+=(double)(clock() - clock_beg) / CLOCKS_PER_SEC;
    //BEFORE
    Vec3_t du;
    
    clock_beg = clock();   
    #pragma omp parallel for schedule (static) private(du) num_threads(Nproc)
    for (int i=0; i<Particles.Size(); i++){
      Particles[i]->x_prev = Particles[i]->x;
      du = (Particles[i]->v + Particles[i]->VXSPH)*deltat*factor;
      Particles[i]->Displacement += du;
      Particles[i]->x += du;
    }

    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (int i=0; i<Particles.Size(); i++){
      Particles[i]->v += Particles[i]->a*deltat/2.*factor;
      //Particles[i]->LimitVel();
    }
    MoveGhost();
    GeneralAfter(*this);
    mov_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;  
    
		clock_beg = clock();
    CalcRateTensors();  //With v and xn+1
    //#ifdef NONLOCK_SUM
    if (nonlock_sum) RateTensorsReduction();
    //#endif
    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (int i=0; i<Particles.Size(); i++){
      //Particles[i]->Mat2Leapfrog(deltat); //Uses density  
      Particles[i]->CalcStressStrain(deltat); //Uses density  
    } 
    stress_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;

		CalcPlasticWorkHeat(deltat);   //Before Thermal increment because it is used
    
    clock_beg = clock();   
    ThermalCalcs(deltat);
    
    thermal_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;

		clock_beg = clock();        
    CalcKinEnergyEqn();    
    CalcIntEnergyEqn();    
    energy_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
    
		steps++;
    if (ct == 30) ct = 0; else ct++;
    


		Time += deltat;
    clock_beg = clock();    
    
    if (isfirst) isfirst = false;    
    
    if (contact){
 		//cout << "checking contact"<<endl;
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

		//if (Time>=tout){
    if (Time>=tout || (double)((clock() - last_output_time) / CLOCKS_PER_SEC) > 60.0) {
      last_output_time = clock();  
      if (Time>=tout ){      
        if (TheFileKey!=NULL) {
          String fn;
          fn.Printf    ("%s_%04d", TheFileKey, idx_out);
          WriteXDMF    (fn.CStr());
          //fn.Printf    ("%s_%.5f", TheFileKey, Time);
          WriteCSV    (fn.CStr());

        }
        idx_out++;
        tout += dtOut;
      }
      
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
			cout << "Max Displacements: "<<max_disp<<endl;
      if (contact) cout<<"Contact Force Sum "<<contact_force_sum<<endl;
      cout << "Int Energy: " << int_energy_sum << ", Kin Energy: " << kin_energy_sum<<endl;
      cout << endl<<"Elapsed since last output: "<< (clock() - last_output_time) / CLOCKS_PER_SEC<<endl;
      
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