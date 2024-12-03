namespace SPH {
////// ACCORDING TO ROBUST AND EFFICIENT MESHFREE SOLID THERMO-MECHANICS
/////////////////// SIMULATION OF FRICTION STIR WELDING  
///////////////////////////////////////
//This calculates all in same time step
//IMPORTANT: x is updated with previous accel!
// 1. x(t+dt) = x +v(t) dt + 1/2 a(t) dt^2
// 2. Calc a(t+dt) using x(t+dt)
// 3. v(t+dt) = v(t) dt + 1/2 (a(t) +a(t+dt)) * dt
///////////////////////////////////////

// // // CalcDensInc(); //TODO: USE SAME KERNEL?
// // // CalcRateTensors();  //With v and xn+1
  // // // Particles[i]->CalcStressStrain(deltat); //Uses density  

// // // CalcAccel(); //Nor density or neither strain rates
  // // // x+= (Particles[i]->v + Particles[i]->VXSPH)*deltat + 0.5 * Particles[i]->a *deltat*deltat;
  // // // Particles[i]->v += Particles[i]->a * deltat;


inline void Domain::SolveDiffUpdateFraser (double tf, double dt, double dtOut, char const * TheFileKey, size_t maxidx) {
	std::cout << "\n--------------Solving---------------------------------------------------------------" << std::endl;
  ostringstream oss_out;

	size_t idx_out = 1;
	double tout = Time;

	//Initializing adaptive time step variables
	deltat = deltatint = deltatmin	= dt;

	clock_t clock_beg;
  double clock_time_spent,start_acc_time_spent, nb_time_spent ,pr_acc_time_spent,acc_time_spent, 
          contact_time_spent, trimesh_time_spent, bc_time_spent,
          mov_time_spent,stress_time_spent,energy_time_spent, dens_time_spent, thermal_time_spent,
          ini_time_spent;
  
  double last_output_time;
  clock_time_spent = contact_time_spent = acc_time_spent = stress_time_spent = energy_time_spent = dens_time_spent = mov_time_spent = thermal_time_spent =
  ini_time_spent = 0.;	

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
  //of << "Displacement, pl_strain, eff_strain_rate, sigma_eq, sigmay, contforcesum"<<endl;
  of << "Time, ";
  for (int cf=0;cf<m_contact_force.size();cf++){
    int sid = contact_surf_id[cf];
    if (cf>0) of <<", ";
      of << "cf"<<sid<<"x, "<<" cf"<<sid<<"y, "<<" cf"<<sid<<"z";
  }
  of <<endl;
  bool check_nb_every_time = false;

  cout << "Main Loop"<<endl;
  cout << "Solver: Fraser "<<endl;
  
  if (nonlock_sum && gradKernelCorr){
    cout << "WARNING: Nishimura summation is not working with Gradient Kernel Correction. Summation changed to Locking." <<endl;
    nonlock_sum = false;
  }
    
  int ct=30;
  std::chrono::duration<double> total_time, elapsed_output;
	auto start_whole = std::chrono::steady_clock::now();  
  
  std::vector <Vec3_t> prev_acc(Particles.Size()); 
  #pragma omp parallel for schedule (static) num_threads(Nproc)	//LUCIANO//LIKE IN DOMAIN->MOVE
  for (int i=0; i<Particles.Size(); i++){
    prev_acc[i] = Particles[i]->a;
  }
  
  cout << std::setprecision(3)<< "Total allocated memory: " <<sizeof(Particle) * Particles.Size() * 1.0e-6 << " MB. "<<endl;
  
  // if (gradKernelCorr){
    // CalcGradCorrMatrix();	}
  last_output_time = clock();   

  while (Time<=tf && idx_out<=maxidx) {
    
    clock_beg = clock();
  
		StartAcceleration(0.);

		double max = 0;
		int imax;
		//#pragma omp parallel for schedule (static) num_threads(Nproc)	//LUCIANO//LIKE IN DOMAIN->MOVE
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
    
		if (!model_damage) {
		if (max > MIN_PS_FOR_NBSEARCH && !isyielding){ //First time yielding, data has not been cleared from first search
			ClearNbData(); 
			MainNeighbourSearch/*_Ext*/();
      CalcPairPosList();  
      SaveNeighbourData();
			if (contact) ContactNbUpdate(this);
			isyielding  = true ;
		}
		}
    
    ini_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;

		if (model_damage && !isfirst) ts_nb_inc = 1; //NEVER SEARCH NBs
		
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

    //after if first
    if (model_damage) {
      CalculateSurface(1);
    }
    
    if (h_update){                
      UpdateSmoothingLength();          
    }
		// //NEW, gradient correction
			if (isfirst) {
				if (gradKernelCorr){
          CalcGradCorrMatrix();	}
				
			}		
		//std::cout << "neighbour_time (chrono, clock): " << clock_time_spent << ", " << neighbour_time.count()<<std::endl;
    
    //cout << "Primary comp "<<endl;
		GeneralBefore(*this);
		PrimaryComputeAcceleration();
    //cout << "done "<<endl;
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
      Particles[i]->Density += deltat*Particles[i]->dDensity;
    } 
    if (dom_bid_type == AxiSymmetric){
      for (int i=0; i<Particles.Size(); i++)
        Particles[i]->Density_real = Particles[i]->Density/(2.0*M_PI*Particles[i]->x(0)); 
    }
    dens_time_spent+=(double)(clock() - clock_beg) / CLOCKS_PER_SEC;

    //cout << "rate tensor"<<endl;
		clock_beg = clock();
    CalcRateTensors();  //With v and xn+1

    //#ifdef NONLOCK_SUM
    if (nonlock_sum)
    RateTensorsReduction();
    //#endif
    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (int i=0; i<Particles.Size(); i++){
      //Particles[i]->Mat2Leapfrog(deltat); //Uses density  
      Particles[i]->CalcStressStrain(deltat); //Uses density  
    } 
    stress_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;		
    
    
    ////// 9. Update Density and Stress: already done
   
    //////////// 11 to 13.
    clock_beg = clock();
    //cout << "Particle 0 accel " << Particles[0]->a<<endl;
    CalcAccel(); //Nor density or neither strain rates
    

    // if (dom_bid_type == AxiSymm_3D)
      // ApplyAxiSymmBC();
    
    //CalcAccelPP();
    //cout << "part 2000 acc "<<Particles[2000]->a<<endl;
    //#ifdef NONLOCK_SUM
    if (nonlock_sum) AccelReduction();
    else {  
      if (dom_bid_type == AxiSymmetric){
      for (int i=0; i<solid_part_count;i++){
        Particles[i]->a *= 2.0 * M_PI;
        Particles[i]->a[0] -= Particles[i]->Sigma(2,2); //WANG Eqn. 40
      }

      cout << "ERROR. Locking Reduction not allowed for Axisymm 2D"<<endl;      
      }
    }
    
    //#endif
		acc_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;
    
    ///// 8. CONTACT FORCES
    cout << "Calc CONTACT"<<endl; 
    clock_beg = clock(); 
    if (contact) {
      if      (contact_alg==Fraser)   CalcContactForces();
      if      (contact_alg==Wang)     CalcContactForcesWang();
      else if (contact_alg==Seo )     CalcContactForces2();
     // else if (contact_alg==LSDyna )  CalcContactForcesLS();
    }
    contact_time_spent +=(double)(clock() - clock_beg) / CLOCKS_PER_SEC;
    cout << "done "<<endl;
    //if (contact) CalcContactForces2();

    // if (contact && !isfirst)
    // #pragma omp parallel for schedule (static) num_threads(Nproc)
    // for (int i=0; i<Particles.Size(); i++)
      // Particles[i]->a += Particles[i]->contforce /*+ Particles[i]->tgforce*/;
    
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
    //CorrectVelAcc(); //CYLINDRICAL SLICE, MOVED TO MOVEGHOST
    //cout << "moving ghost "<<endl;
    MoveGhost(); 
    //cout << " ghost doe "<<endl;    
   #pragma omp parallel for schedule (static) private(du) num_threads(Nproc)
    for (int i=0; i<Particles.Size(); i++){
      Particles[i]->x_prev = Particles[i]->x;
      du = (Particles[i]->v + Particles[i]->VXSPH)*deltat + 0.5 * Particles[i]->a *deltat*deltat;
      Particles[i]->Displacement += du;
      Particles[i]->x += du;
    }
    mov_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;  
    
    
		
    double factor = 1.;

    clock_beg = clock();

    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (int i=0; i<Particles.Size(); i++){
      Particles[i]->v += Particles[i]->a * deltat;
    }

    //CorrectVelAcc(); //CYLINDRICAL SLICE, MOVED TO MOVEGHOST
    MoveGhost();   

    GeneralAfter(*this);//Reinforce BC vel   
    mov_time_spent += (double)(clock() - clock_beg) / CLOCKS_PER_SEC;  

    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (int i=0; i<Particles.Size(); i++)
      prev_acc[i] = Particles[i]->a;

		if (model_damage) CalcDamage();
		
		CalcPlasticWorkHeat(deltat);   //Before Thermal increment because it is used
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
    
    //auto last_output_time = start_whole;
    
    
    cout << "time "<< (double)((clock() - last_output_time) / CLOCKS_PER_SEC)<<endl;
		//if (Time>=tout || std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - last_output_time).count() > 60.0){
    //cout << "elapsed "<<(double)((clock() - last_output_time) / CLOCKS_PER_SEC)<<endl;
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
      oss_out.str("");
			oss_out << "\n---------------------------------------\n Total CPU time: "<<total_time.count() << endl;
      double acc_time_spent_perc = acc_time_spent/total_time.count();

      oss_out << std::setprecision(2);
      oss_out << "step number: "<<steps <<endl;
      oss_out << "Calculation Times\nAccel: "<<acc_time_spent_perc<<"%,  ";
      oss_out << "Density: "<<dens_time_spent/total_time.count()<<"%,  ";
      oss_out << "Stress: "  <<stress_time_spent/total_time.count()<<"%,  "<<endl;
      oss_out << "Energy: "  <<energy_time_spent/total_time.count()<<"%,  ";
      oss_out << "Contact: " <<contact_time_spent/total_time.count()<<"%,  ";
      oss_out << "Nb: "      <<nb_time_spent/total_time.count()<<"%,  ";
      oss_out << "Update: " <<mov_time_spent/total_time.count()<<"%,  ";
      oss_out <<endl;
      
			oss_out << "Output No. " << idx_out << " at " << Time << " has been generated" << std::endl;
			oss_out << "Current Time Step = " <<deltat<<std::endl;
			oss_out << "Max plastic strain: " <<max<< "in particle" << imax << endl;
      if (max > 0.)
        oss_out<<"Plastic Work "<<plastic_work<<endl;
      oss_out.precision(6);
			oss_out << "Max Displacements (No Cont Surf): "<<max_disp<<endl;
      if (contact) {
        oss_out<<"Contact Force Sum "<<contact_force_sum<<", Reaction Sum "<< contact_reaction_sum<<endl;
        oss_out<<"Contact Friction Work "<<contact_friction_work<<endl;
        oss_out<<"External Forces Work "<< ext_forces_work<<endl;
        if (cont_heat_cond)
          oss_out << "Total contact heat flux" << accum_cont_heat_cond <<endl;
      }
      oss_out << "Int Energy: " << int_energy_sum << ", Kin Energy: " << kin_energy_sum<<endl;
      if (thermal_solver)
        oss_out << "Max, Min, Avg temps: "<< m_maxT << ", " << m_minT << ", " << (m_maxT+m_minT)/2. <<std::endl;    

      ofprop <<getTime() << ", "<<m_scalar_prop<<endl;
			
      // for (int p=0;p<Particles.Size();p++){
				// if (Particles[p]->print_history)
          // of << Particles[p]->Displacement << ", "<<Particles[p]->pl_strain<<", "<<Particles[p]->eff_strain_rate<<", "<< 
          // Particles[p]->Sigma_eq<<", "  <<  Particles[p]->Sigmay << ", " <<
          // contact_force_sum << endl;
			// }
      of <<Time<<", ";
      for (int cf=0;cf<m_contact_force.size();cf++){
        if (cf>0)of <<", ";
        of << m_contact_force[cf](0)<<", "<<m_contact_force[cf](1)<<", "<<m_contact_force[cf](2);
      }
      of <<endl;
			if (model_damage){
				double max_dam = 0.;
				int dam_count =0;
				for (int p=0;p<Particles.Size();p++){
					if (Particles[p]->dam_D>max_dam) 
						max_dam = Particles[p]->dam_D;
					if (Particles[p]->dam_D==1.0) 
						dam_count++;
				}
			cout << "Max damage is "<<max_dam<<", full damage count: "<<dam_count<<endl;
			}
      
      oss_out << "-----------------------------------------------------------------------------"<<endl;
      cout << oss_out.str();
      out_file << oss_out.str();
      
		}
    // if (auto_ts)      CheckMinTSVel();
    // if (auto_ts_acc)  CheckMinTSAccel();
    // if (auto_ts || auto_ts_acc)  AdaptiveTimeStep();
    
	
	}


	of.close(); //History 
  ofprop.close(); //Scalar prop
	
	std::cout << "\n--------------Solving is finished---------------------------------------------------" << std::endl;

}

};
