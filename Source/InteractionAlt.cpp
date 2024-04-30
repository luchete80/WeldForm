#include "Domain.h"

///////////////////////////////////////////////////////////////////
/*  NEW FUNCTION TO CALCULATE ACCELERATION, REMOVING STRAIN AND ROTATION RATES
*** AND ALSO DENSITY; IN ORDER TO CALCULATE THEM SEPARATELY *//////
//  NOTE: ONLY FOR FREE PARTICLES

#define ID_TEST 2000
namespace SPH{
  
  
inline void Domain::CalcAccel() {
  Particle *P1, *P2;
  double dam_f; //if not damage
	
  #pragma omp parallel for schedule (static) private (P1,P2,dam_f) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t k=0; k<Nproc;k++) 
	#else
	for (int k=0; k<Nproc;k++) 
	#endif	
	{
  for (size_t p=0; p<SMPairs[k].Size();p++) {
    
    //#ifndef NONLOCK_SUM
    if (!nonlock_sum){
    P1	= Particles[SMPairs[k][p].first];
    P2	= Particles[SMPairs[k][p].second];	
    //#else
    } else {
    P1 = Particles[std::min(SMPairs[k][p].first, SMPairs[k][p].second)];
    P2 = Particles[std::max(SMPairs[k][p].first, SMPairs[k][p].second)];
    //#endif
    }
    double h	= (P1->h+P2->h)/2;
    Vec3_t xij	= P1->x - P2->x;
    dam_f = 1.0;
    //Periodic_X_Correction(xij, h, P1, P2);

    double rij	= norm(xij);
    
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
		double Alpha	= (P1->Alpha + P2->Alpha)/2.0;
		double Beta	= (P1->Beta + P2->Beta)/2.0;

    //In axi-Symm mode, real density is only considered
    
    di = P1->Density;
    mi = P1->Mass;

    dj = P2->Density;
    mj = P2->Mass;
    
    if (model_damage) {
      if (dam_D[k][p] > 0.0 ) dam_f = 1.0 - dam_D[k][p];
      // if (dam_D[k][p]     
    }
		
    if (dam_f > 0.0) {
      if (dam_f > 1.0 ) dam_f = 1.0;
		Vec3_t vij	= P1->v - P2->v;
		double GK	= GradKernel(Dimension, KernelType, rij/h, h);
		double K	= Kernel(Dimension, KernelType, rij/h, h);
		
		// double GK	= m_kernel.gradW(rij/h);
		// double K		= m_kernel.W(rij/h);
		
    //Real density for sound speed calc

    //m_clock_begin = clock();
		// Artificial Viscosity
    Mat3_t PIij;
    set_to_zero(PIij);
    if (Alpha!=0.0 || Beta!=0.0)
    {
      double MUij = h*dot(vij,xij)/(rij*rij+0.01*h*h);					///<(2.75) Li, Liu Book
      double Cij;
      double Ci,Cj;
      if (dom_bid_type == AxiSymmetric){ //CALCULATED DENSITY
        di/=(2.0*M_PI*P1->x(0));
        dj/=(2.0*M_PI*P2->x(0));
      }
      if (!P1->IsFree) Ci = SoundSpeed(P2->PresEq, P2->Cs, di, P2->RefDensity); 
      else            Ci = SoundSpeed(P1->PresEq, P1->Cs, di, P1->RefDensity);
      if (!P2->IsFree) Cj = SoundSpeed(P1->PresEq, P1->Cs, dj, P1->RefDensity); else Cj = SoundSpeed(P2->PresEq, P2->Cs, dj, P2->RefDensity);
      Cij = 0.5*(Ci+Cj);
      ////Mod density for artif visc calc
      if (dom_bid_type == AxiSymmetric){
        di*=(2.0*M_PI*P1->x(0));
        dj*=(2.0*M_PI*P2->x(0));
      }      
      if (dot(vij,xij)<0) PIij = (Alpha*Cij*MUij+Beta*MUij*MUij)/(0.5*(di+dj)) * I;		///<(2.74) Li, Liu Book
      //cout << "PIij"<<PIij<<endl;
    }
    //m_forces_artifvisc_time += (double)(clock() - m_clock_begin) / CLOCKS_PER_SEC;
    
    Mat3_t Sigmaj,Sigmai;
    set_to_zero(Sigmaj);
    set_to_zero(Sigmai);
    Sigmai = P1->Sigma;
    Sigmaj = P2->Sigma;

    // Tensile Instability
    Mat3_t TIij;
    set_to_zero(TIij);
    if (P1->TI > 0.0 || P2->TI > 0.0) {
      //cout << "P1->TIR" << P1->TIR<<endl;
      TIij = pow((K/Kernel(Dimension, KernelType, (P1->TIInitDist + P2->TIInitDist)/(2.0*h), h)),(P1->TIn+P2->TIn)/2.0)*(P1->TIR+P2->TIR);
      //TIij = pow((K/m_kernel.W((P1->TIInitDist + P2->TIInitDist)/(2.0*h))),(P1->TIn+P2->TIn)/2.0)*(P1->TIR+P2->TIR);
    }
		Mat3_t StrainRate,RotationRate;
		set_to_zero(StrainRate);
		set_to_zero(RotationRate);
		
		//NEW
		Mat3_t GKc[2];
		double m, mc[2];
		GKc[0] = GK * P1->gradCorrM;
		GKc[1] = GK * P2->gradCorrM;
		if (gradKernelCorr){
		}
		// Calculating the forces for the particle 1 & 2
		Vec3_t temp = 0.0;
		double temp1 = 0.0;
		Vec3_t temp_c[2];
		double temp1_c[2];
		Vec3_t vc[2];
		
		if (gradKernelCorr){
			for (int i=0;i<2;i++){
				Mult (GKc[i], xij, vc[i]);
			}
		}

		// NEW
		if (!gradKernelCorr) {
      if (dom_bid_type != AxiSymmetric){
        if (GradientType == 0)
          Mult( GK*xij , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij) , temp);
        else if (GradientType == 1)
          Mult( GK*xij , ( 1.0/(di*dj)*(Sigmai + Sigmaj)           + PIij + TIij) , temp);
        else 
          Mult( GK*xij , ( 1.0/(di*dj)*(Sigmai - Sigmaj)           + PIij + TIij) , temp); ////////// SEEE CAMPBELL 2000
      } else { //Only w/o gradient corr
        //if (GradientType == 0 ){
          //etaDens
          //Wang Eqn 40 and Joshi Eqn 27
          Vec3_t wij = GK*xij;
          // di = P1->etaDens;
          // dj = P2->etaDens;
          Vec3_t av;
          Mult (wij, PIij,av);
          //  cout << "AVISC TERM "<<av<<endl;
          temp[0] = (Sigmai(0,0)*P1->x(0)/(di*di) + Sigmaj(0,0) *P2->x(0)/(dj*dj)) *wij(0) + 
                    (Sigmai(0,1)*P1->x(0)/(di*di) + Sigmaj(0,1) *P2->x(0)/(dj*dj)) *wij(1) +PIij;  ////dvr/dt 2PI can go in the reduction

          temp[1] = (Sigmai(0,1)*P1->x(0)/(di*di) + Sigmaj(0,1) *P2->x(0)/(dj*dj)) *wij(0) + 
                    (Sigmai(1,1)*P1->x(0)/(di*di) + Sigmaj(1,1) *P2->x(0)/(dj*dj)) *wij(1) ;  ////dvr/dt 2PI can go in the reduction
          
          
        //}
      }//axisymm
		} else {
				//Should be replaced  dot( xij , GK*xij ) by dot( xij , v )
				//Left in vector form and multiply after??
        if (GradientType == 0){
          for (int i=0;i<2;i++){
            Mult( vc[i] , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , temp_c[i]);
          }
        } else if (GradientType == 1) {
          for (int i=0;i<2;i++){
            Mult( vc[i] , ( 1.0/(di*dj)*(Sigmai + Sigmaj)           + PIij + TIij ) , temp_c[i]);
            } 
          } else if (GradientType == 2) {
          for (int i=0;i<2;i++){
            Mult( vc[i] , ( 1.0/(di*dj)*(Sigmai - Sigmaj)           + PIij + TIij ) , temp_c[i]);
          }          
        }
		}//Grad Corr
    

		if (Dimension == 2 /*&& dom_bid_type != AxiSymmetric*/) temp(2) = 0.0; //PLANE STRAIN
		
		if (!gradKernelCorr){
			temp1 = dot( vij , GK*xij );
		} else {
			for (int i=0;i<2;i++){			//TODO: DO THIS ONCE!
				temp1_c[i] = dot( vij , vc[i] );
			}
		}
    
    //#ifdef NONLOCK_SUM
    if (nonlock_sum)
    //if (!gradKernelCorr) 
    pair_force[first_pair_perproc[k] + p] = temp; //SHOULD ALSO MULTIPLY ACCEL AFTER
    
    // if (SMPairs[k][p].first == ID_TEST || SMPairs[k][p].second == ID_TEST)
      // cout << "i j temp mj: "<<SMPairs[k][p].first<<", "<<SMPairs[k][p].second<<", "<< temp;
    // if (SMPairs[k][p].first == ID_TEST) cout << mj <<", ";
    // else if (SMPairs[k][p].second == ID_TEST) cout << mi<<", ";
    
    // if (SMPairs[k][p].first == ID_TEST) cout << "-"<<endl;
    // else if (SMPairs[k][p].second == ID_TEST) cout << "+" <<endl;
    //#else  ////NONLOCK
    else { 
		// Locking the particle 1 for updating the properties
		omp_set_lock(&P1->my_lock);
			if (!gradKernelCorr){
				P1->a					+= dam_f*mj * temp;
				//P1->dDensity	+= mj * (di/dj) * temp1;
			} else{
				P1->a					+= dam_f *mj * temp_c[0];
				//P1->dDensity	+= mj * (di/dj) * temp1_c[0];
			}

		omp_unset_lock(&P1->my_lock);

		// Locking the particle 2 for updating the properties
		omp_set_lock(&P2->my_lock);
			if (!gradKernelCorr){
				P2->a					-= dam_f * mi * temp;				
			}else {
				P2->a					-= dam_f * mi * temp_c[1];
			}
		omp_unset_lock(&P2->my_lock);
    //#endif
    }//nonlock_sum
    
    }//dam_f
  }//MAIN FOR IN PAIR
  }//MAIN FOR PROC

}

inline void Domain::AccelReduction(){
  if (solid_part_count > 0){
    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (int i=0; i<solid_part_count;i++)
      Particles[i]->a = 0.;
    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (int i=0; i<solid_part_count;i++){
      for (int n=0;n<ipair_SM[i];n++){  
        Particles[i]->a += Particles[Anei[i][n]]->Mass * pair_force[Aref[i][n]];}
      for (int n=0;n<jpair_SM[i];n++){   
        Particles[i]->a -= Particles[Anei[i][MAX_NB_PER_PART-1-n]]->Mass * pair_force[Aref[i][MAX_NB_PER_PART-1-n]];}    
    }
    if (dom_bid_type == AxiSymmetric){
      //ADD HOOP ACCEL AND MULT BY 2PI
      #pragma omp parallel for schedule (static) num_threads(Nproc)
      for (int i=0; i<solid_part_count;i++){
        Particles[i]->a *= 2.0 * M_PI; //PREVIOUSLY
        Particles[i]->a[0] -= 2.0 * M_PI * Particles[i]->Sigma(2,2) / Particles[i]->Density; //WANG Eqn. 40
      }
    }
  } else {
    cout << "ERROR, particle count not defined"<<endl;
  }
}

//Similar but not densities
inline void Domain::CalcRateTensors() {
  Particle *P1, *P2;          
	#pragma omp parallel for schedule (static) private (P1,P2) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t k=0; k<Nproc;k++) 
	#else
	for (int k=0; k<Nproc;k++) 
	#endif	
	{
  for (size_t p=0; p<SMPairs[k].Size();p++) {
    //#ifndef NONLOCK_SUM
    if (!nonlock_sum){
    P1	= Particles[SMPairs[k][p].first];
    P2	= Particles[SMPairs[k][p].second];	
    } else {
    //#else
    P1 = Particles[std::min(SMPairs[k][p].first, SMPairs[k][p].second)];
    P2 = Particles[std::max(SMPairs[k][p].first, SMPairs[k][p].second)];
    //#endif
    }
    double h	= (P1->h+P2->h)/2;
    Vec3_t xij	= P1->x - P2->x;

	//Periodic_X_Correction(xij, h, P1, P2);

	double rij	= norm(xij);

  double clock_begin;

	// if ((rij/h)<=Cellfac)
	// {
  double di=0.0,dj=0.0,mi=0.0,mj=0.0;

    di = P1->Density;
    mi = P1->Mass;

    dj = P2->Density;
    mj = P2->Mass;

		Vec3_t vij	= P1->v - P2->v;
    
    int i1 = std::min(SMPairs[k][p].first, SMPairs[k][p].second);
		int i2 = std::max(SMPairs[k][p].first, SMPairs[k][p].second);
    
		double GK	= GradKernel(Dimension, KernelType, rij/h, h);
		double K	= Kernel(Dimension, KernelType, rij/h, h);

		Mat3_t Sigmaj,Sigmai;
		set_to_zero(Sigmaj);
		set_to_zero(Sigmai);
		Sigmai = P1->Sigma;
		Sigmaj = P2->Sigma;

		// NoSlip BC velocity correction
		Vec3_t vab = vij;

		Mat3_t StrainRate,RotationRate;
		set_to_zero(StrainRate);
		set_to_zero(RotationRate);
		
		//NEW
		Mat3_t GKc[2];
		double m, mc[2];
		GKc[0] = GK * P1->gradCorrM;
		GKc[1] = GK * P2->gradCorrM;
		if (gradKernelCorr){
		}

    //m_clock_begin = clock();		
		
    Mat3_t StrainRate_c[2],RotationRate_c[2]; //Corrected gradients

		// // // // Calculation strain rate tensor
		// // // //ORIGINAL FORM			
		//if (!gradKernelCorr){
    StrainRate(0,0) = 2.0*vab(0)*xij(0);
    StrainRate(0,1) = vab(0)*xij(1)+vab(1)*xij(0);
    StrainRate(0,2) = vab(0)*xij(2)+vab(2)*xij(0);
    StrainRate(1,0) = StrainRate(0,1);
    StrainRate(1,1) = 2.0*vab(1)*xij(1);
    StrainRate(1,2) = vab(1)*xij(2)+vab(2)*xij(1);
    // if (dom_bid_type == AxiSymmetric){
      // StrainRate(0,2) = StrainRate(1,2) = 0.;
      // RotationRate(0,2) = RotationRate(1,2) = 0.;
    // }
    StrainRate(2,0) = StrainRate(0,2);
    StrainRate(2,1) = StrainRate(1,2);
    StrainRate(2,2) = 2.0*vab(2)*xij(2);
    StrainRate	= -0.5 * GK * StrainRate;
    

    
    RotationRate(0,1) = vab(0)*xij(1)-vab(1)*xij(0);
    RotationRate(0,2) = vab(0)*xij(2)-vab(2)*xij(0);
    RotationRate(1,2) = vab(1)*xij(2)-vab(2)*xij(1);
    RotationRate(1,0) = -RotationRate(0,1);
    RotationRate(2,0) = -RotationRate(0,2);
    RotationRate(2,1) = -RotationRate(1,2);
    RotationRate	  = -0.5 * GK * RotationRate;
    
			// if (StrainRate(2,2)<-1.e-3)
        
			Mat3_t gradv[2],gradvT[2];
			
			//cout<<"gradv"<<gradv[0]<<endl;
			Vec3_t gradK; 
			Mult(GK * P1->gradCorrM,xij,gradK);
			Dyad (vab,gradK,gradv[0]); //outer product. L, velocity gradient tensor
			Mult(GK * P2->gradCorrM,xij,gradK);
			Dyad (vab,gradK,gradv[1]); //outer product. L, velocity gradient tensor
			
			for (int j=0;j<2;j++){
				Trans(gradv[j],gradvT[j]);
				StrainRate_c[j] 	= -0.5*(gradv[j] + gradvT[j]);
				RotationRate_c[j] = -0.5*(gradv[j] - gradvT[j]);
			}
      
		// Calculating the forces for the particle 1 & 2
		Vec3_t temp = 0.0;
		double temp1 = 0.0;
		Vec3_t temp_c[2];
		double temp1_c[2];
		Vec3_t vc[2];
		
		if (gradKernelCorr){
			for (int j=0;j<2;j++){
				Mult (GKc[j], xij, vc[j]);
			}
		}


		if (Dimension == 2) temp(2) = 0.0;
		
		if (!gradKernelCorr){
			temp1 = dot( vij , GK*xij );
		} else {
			for (int i=0;i<2;i++){			//TODO: DO THIS ONCE!
				temp1_c[p] = dot( vij , vc[p] );
			}
		}
    
    clock_begin = clock();
		// Locking the particle 1 for updating the properties

    //#ifdef NONLOCK_SUM
    if (nonlock_sum){
    //if (!gradKernelCorr) 
    pair_StrainRate[first_pair_perproc[k] + p] = StrainRate; //SHOULD ALSO MULTIPLY ACCEL AFTER
    pair_RotRate[first_pair_perproc[k] + p] = RotationRate; //SHOULD ALSO MULTIPLY ACCEL AFTER
    //#else
    } else {
		omp_set_lock(&P1->my_lock);

      float mj_dj= mj/dj;

      if (!gradKernelCorr){
        P1->StrainRate 		= P1->StrainRate + mj_dj*StrainRate;
        P1->RotationRate 	= P1->RotationRate + mj_dj*RotationRate;
      }
      else {
        P1->StrainRate 		= P1->StrainRate 		+ mj_dj * StrainRate_c[0];
        P1->RotationRate 	= P1->RotationRate 	+ mj_dj * RotationRate_c[0];
      }

		omp_unset_lock(&P1->my_lock);

		// Locking the particle 2 for updating the properties
		omp_set_lock(&P2->my_lock);
	
      float mi_di = mi/di;
      if (!gradKernelCorr){
        P2->StrainRate	 = P2->StrainRate + mi_di*StrainRate;
        P2->RotationRate = P2->RotationRate + mi_di*RotationRate;
      } else {
        P2->StrainRate = P2->StrainRate 		+ mi_di*StrainRate_c[1];
        P2->RotationRate = P2->RotationRate + mi_di*RotationRate_c[1];
      }

		omp_unset_lock(&P2->my_lock);
    //#endif
    }//nonlock_sum
    }//FOR PAIRS
  }//FOR NPROC
}


// TODO: TEMPLATIZE, at least by type, by Reduction double, 
inline void Domain::RateTensorsReduction(){
  //Not necesay to set to zero here. Are in domain
  #pragma omp parallel for schedule (static) num_threads(Nproc)
  for (int i=0; i<solid_part_count;i++){
    for (int n=0;n<ipair_SM[i];n++){    
      double mjdj = Particles[Anei[i][n]]->Mass /Particles[Anei[i][n]]->Density;
      Particles[i]->StrainRate    = Particles[i]->StrainRate   + mjdj * pair_StrainRate[Aref[i][n]];
      Particles[i]->RotationRate  = Particles[i]->RotationRate + mjdj * pair_RotRate[Aref[i][n]];      
    }
    for (int n=0;n<jpair_SM[i];n++){   
      double mjdj = Particles[Anei[i][MAX_NB_PER_PART-1-n]]->Mass / Particles[Anei[i][MAX_NB_PER_PART-1-n]]->Density;
      Particles[i]->StrainRate    = Particles[i]->StrainRate   + mjdj * pair_StrainRate[Aref[i][MAX_NB_PER_PART-1-n]];
      Particles[i]->RotationRate  = Particles[i]->RotationRate + mjdj * pair_RotRate[Aref[i][MAX_NB_PER_PART-1-n]];
    } 
  }

  if (dom_bid_type == AxiSymmetric){ //WANG 2015 EQN 17, DIRECT HOOP
    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (int i=0; i<solid_part_count;i++){
      double psi;
      //if (Particles[i]->x(0)*Particles[i]->x(0) < Particles[i]->h*Particles[i]->h){
        if (Particles[i]->x(0) > 0.0) psi =  0.01 * Particles[i]->h;
        else                          psi = -0.01 * Particles[i]->h;
        Particles[i]->StrainRate(2,2) = Particles[i]->v(0)/(psi + Particles[i]->x(0)); //DIRECT HOOP STRAIN RATE, Wang 2015
        Particles[i]->StrainRate(0,2) = 0.0; Particles[i]->StrainRate(2,0) = 0.0;
         Particles[i]->RotationRate(0,2) = 0.0; Particles[i]->RotationRate(2,0) = 0.0;
      //}    
    }
  }
}


// TODO: USED CALCULATED KERNELKS
inline void Domain::CalcDensInc() {
	double dam_f;
  Particle *P1, *P2;
	#pragma omp parallel for schedule (static) private (P1,P2,dam_f) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t k=0; k<Nproc;k++) 
	#else
	for (int k=0; k<Nproc;k++) 
	#endif	
	{
    for (size_t p=0; p<SMPairs[k].Size();p++) {
      if (!nonlock_sum){
      //#ifndef NONLOCK_SUM
      P1	= Particles[SMPairs[k][p].first];
      P2	= Particles[SMPairs[k][p].second];	
      //#else
      } else {
      P1 = Particles[std::min(SMPairs[k][p].first, SMPairs[k][p].second)];
      P2 = Particles[std::max(SMPairs[k][p].first, SMPairs[k][p].second)];
      //#endif
      }
      dam_f = 1.0;
      double h	= (P1->h+P2->h)/2;
      Vec3_t xij	= P1->x - P2->x;

      //Periodic_X_Correction(xij, h, P1, P2);

      double rij	= norm(xij);

      // if ((rij/h)<=Cellfac)
      // {
      double di=0.0,dj=0.0,mi=0.0,mj=0.0;

      di = P1->Density;
      mi = P1->Mass;

      dj = P2->Density;
      mj = P2->Mass;

      Vec3_t vij	= P1->v - P2->v;
      
      double GK	= GradKernel(Dimension, KernelType, rij/h, h);
      double K	= Kernel(Dimension, KernelType, rij/h, h);
    
      //NEW
      Mat3_t GKc[2];
      double m, mc[2];
      GKc[0] = GK * P1->gradCorrM;
      GKc[1] = GK * P2->gradCorrM;
      if (gradKernelCorr){
      }
      // Calculating the forces for the particle 1 & 2
      Vec3_t temp = 0.0;
      double temp1 = 0.0;
      Vec3_t temp_c[2];
      double temp1_c[2];
      Vec3_t vc[2];
      
      if (gradKernelCorr){
        for (int i=0;i<2;i++){
          Mult (GKc[i], xij, vc[i]);
        }
      }


      if (Dimension == 2) temp(2) = 0.0;
      
      if (!gradKernelCorr){
        temp1 = dot( vij , GK*xij );
      } else {
        for (int i=0;i<2;i++){			//TODO: DO THIS ONCE!
          temp1_c[i] = dot( vij , vc[i] );
        }
      }
      //#ifdef NONLOCK_SUM
      if (model_damage) {
        if (dam_D[k][p] > 0.0 ) dam_f = 1.0 - dam_D[k][p];
        // if (dam_D[k][p]     
      } 
      if (dam_f > 0.0){
        if (dam_f>1.0) dam_f = 1.0;
      if (nonlock_sum)
        pair_densinc[first_pair_perproc[k] + p] = temp1;
      //#else
      else {
      // Locking the particle 1 for updating the properties
      if (dom_bid_type != AxiSymmetric) {
      omp_set_lock(&P1->my_lock);
        if (!gradKernelCorr){
          P1->dDensity	+= dam_f * mj * (di/dj) * temp1;
        } else{
          P1->dDensity	+= dam_f *mj * (di/dj) * temp1_c[0];
        }		
      omp_unset_lock(&P1->my_lock);

      // Locking the particle 2 for updating the properties
      omp_set_lock(&P2->my_lock);
        if (!gradKernelCorr){
          P2->dDensity	+= dam_f * mi * (dj/di) * temp1;							
        }else {
          P2->dDensity	+= dam_f * mi * (dj/di) * temp1_c[1];
        }
      omp_unset_lock(&P2->my_lock);
      
      } else{ //HERE DENSITY IS NOT STANDARD DENSITY BUT : 2*PI*r*rho WANG EQN 56
      ///// AISYMM NOW WORKING WITH STD REDUCTION (THIS)
        omp_set_lock(&P1->my_lock);
          if (!gradKernelCorr){
            P1->dDensity	+= dam_f * mj * temp1;
          } else{
            P1->dDensity	+= dam_f * mj * temp1_c[0];
          }		
        omp_unset_lock(&P1->my_lock);        

          if (!gradKernelCorr){
            P2->dDensity	+= dam_f * mi * temp1;
          } else{
            P2->dDensity	+= dam_f * mi * temp1_c[1];
          }		
        omp_unset_lock(&P1->my_lock); 
      }

      //#endif
      }//nonlock_sum
      } //if dam_f > 0.0
    }//FOR PAIRS
  }//FOR NPROC

}

inline void Domain::DensReduction(){
  if (dom_bid_type != AxiSymmetric) {
    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (int i=0; i<solid_part_count;i++){
      Particles[i]->dDensity = 0.;
      for (int n=0;n<ipair_SM[i];n++){ 
        Particles[i]->dDensity += Particles[Anei[i][n]]->Mass /Particles[Anei[i][n]]->Density * pair_densinc[Aref[i][n]];
      }
      for (int n=0;n<jpair_SM[i];n++){   
        Particles[i]->dDensity += Particles[Anei[i][MAX_NB_PER_PART-1-n]]->Mass / Particles[Anei[i][MAX_NB_PER_PART-1-n]]->Density * pair_densinc[Aref[i][MAX_NB_PER_PART-1-n]];    
      }      
      Particles[i]->dDensity *= Particles[i]->Density;
    }
  } else { //HERE DENSITY IS NOT STANDARD DENSITY BUT : 2*PI*r*rho WANG EQN 56
    #pragma omp parallel for schedule (static) num_threads(Nproc)
    for (int i=0; i<solid_part_count;i++){
      Particles[i]->dDensity = 0.;
      for (int n=0;n<ipair_SM[i];n++){ 
        Particles[i]->dDensity += Particles[Anei[i][n]]->Mass * pair_densinc[Aref[i][n]];
      }
      for (int n=0;n<jpair_SM[i];n++){   
        Particles[i]->dDensity += Particles[Anei[i][MAX_NB_PER_PART-1-n]]->Mass * pair_densinc[Aref[i][MAX_NB_PER_PART-1-n]];    
      }    
      //Particles[i]->dDensity *= Particles[i]->Density;
    }    
  }
  
}

inline void Domain::CalcForceSOA(int &i,int &j) {

  Particle * P1, *P2;
	double h	= (P1->h+P2->h)/2;
	Vec3_t xij	= P1->x - P2->x;

	Periodic_X_Correction(xij, h, P1, P2);

	double rij	= norm(xij);

  double clock_begin;

	// if ((rij/h)<=Cellfac)
	// {
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
		double Alpha	= (P1->Alpha + P2->Alpha)/2.0;
		double Beta	= (P1->Beta + P2->Beta)/2.0;

		// if (!P1->IsFree) {
			// di = DensitySolid(P2->PresEq, P2->Cs, P2->P0,P1->Pressure, P2->RefDensity);
			// mi = P1->FPMassC * P2->Mass;
		// } else {
			// di = P1->Density;
			// mi = P1->Mass;
		// }
		// if (!P2->IsFree) {
			// dj = DensitySolid(P1->PresEq, P1->Cs, P1->P0,P2->Pressure, P1->RefDensity);
			// mj = P2->FPMassC * P1->Mass;
		// } else {
			// dj = P2->Density;
			// mj = P2->Mass;
		// }
		

		Vec3_t vij	= P1->v - P2->v;
		
		double GK	= GradKernel(Dimension, KernelType, rij/h, h);
		double K	= Kernel(Dimension, KernelType, rij/h, h); 
    
		// Artificial Viscosity
		Mat3_t PIij;
		set_to_zero(PIij);
		if (Alpha!=0.0 || Beta!=0.0)
		{
			double MUij = h*dot(vij,xij)/(rij*rij+0.01*h*h);					///<(2.75) Li, Liu Book
			double Cij;
			double Ci,Cj;
			if (!P1->IsFree)  Ci = SoundSpeed(P2->PresEq, P2->Cs, di, P2->RefDensity); 
      else              Ci = SoundSpeed(P1->PresEq, P1->Cs, di, P1->RefDensity);
			if (!P2->IsFree)  Cj = SoundSpeed(P1->PresEq, P1->Cs, dj, P1->RefDensity); 
      else              Cj = SoundSpeed(P2->PresEq, P2->Cs, dj, P2->RefDensity);
			Cij = 0.5*(Ci+Cj);
			
			if (dot(vij,xij)<0) PIij = (Alpha*Cij*MUij+Beta*MUij*MUij)/(0.5*(di+dj)) * I;		///<(2.74) Li, Liu Book
		}
    //m_forces_artifvisc_time += (double)(clock() - m_clock_begin) / CLOCKS_PER_SEC;

		Mat3_t Sigmaj,Sigmai;
		set_to_zero(Sigmaj);
		set_to_zero(Sigmai);
		//Sigmai = P1->Sigma;
    //Sigmaj = P2->Sigma;
    
    double tempi[6],tempj[6];
		for (int k=0;k<6;k++){ //First the diagonal
			tempi[k]=sigma[6*i+k];
			tempj[k]=sigma[6*j+k];
		}
		
		Sigmai = FromFlatSym(tempi);
		Sigmaj = FromFlatSym(tempj);
    
		

//		if (P1->IsFree) Sigmai = P1->Sigma; else  Sigmai = P2->Sigma;
//		if (P2->IsFree) Sigmaj = P2->Sigma; else  Sigmaj = P1->Sigma;

		// Tensile Instability
		Mat3_t TIij;
		set_to_zero(TIij);
		if (P1->TI > 0.0 || P2->TI > 0.0) 
			TIij = pow((K/Kernel(Dimension, KernelType, (P1->TIInitDist + P2->TIInitDist)/(2.0*h), h)),(P1->TIn+P2->TIn)/2.0)*(P1->TIR+P2->TIR);
			//TIij = pow((K/m_kernel.W((P1->TIInitDist + P2->TIInitDist)/(2.0*h))),(P1->TIn+P2->TIn)/2.0)*(P1->TIR+P2->TIR);

		// NoSlip BC velocity correction
		Vec3_t vab = 0.0;
		if (P1->IsFree*P2->IsFree) {
			vab = vij;
		} else {
			if (P1->NoSlip || P2->NoSlip) {
				// No-Slip velocity correction
				if (P1->IsFree)	vab = P1->v - (2.0*P2->v-P2->NSv); else vab = (2.0*P1->v-P1->NSv) - P2->v;
			}
			// Please check
			if (!(P1->NoSlip || P2->NoSlip)) {
				if (P1->IsFree) vab = P1->v - P2->vb; else vab = P1->vb - P2->v;
//				if (P1->IsFree) vab(0) = P1->v(0) + P2->vb(0); else vab(0) = -P1->vb(0) - P2->v(0);
			}
		}

		Mat3_t StrainRate,RotationRate;
		set_to_zero(StrainRate);
		set_to_zero(RotationRate);
		
		//NEW
		Mat3_t GKc[2];
		double m, mc[2];
		GKc[0] = GK * P1->gradCorrM;
		GKc[1] = GK * P2->gradCorrM;
		if (gradKernelCorr){
		}

    //m_clock_begin = clock();		
		
    Mat3_t StrainRate_c[2],RotationRate_c[2]; //Corrected gradients

		// // // // Calculation strain rate tensor
		// // // //ORIGINAL FORM			
		//if (!gradKernelCorr){
			StrainRate(0,0) = 2.0*vab(0)*xij(0);
			StrainRate(0,1) = vab(0)*xij(1)+vab(1)*xij(0);
			StrainRate(0,2) = vab(0)*xij(2)+vab(2)*xij(0);
			StrainRate(1,0) = StrainRate(0,1);
			StrainRate(1,1) = 2.0*vab(1)*xij(1);
			StrainRate(1,2) = vab(1)*xij(2)+vab(2)*xij(1);
			StrainRate(2,0) = StrainRate(0,2);
			StrainRate(2,1) = StrainRate(1,2);
			StrainRate(2,2) = 2.0*vab(2)*xij(2);
			StrainRate	= -0.5 * GK * StrainRate;
      
			
			RotationRate(0,1) = vab(0)*xij(1)-vab(1)*xij(0);
			RotationRate(0,2) = vab(0)*xij(2)-vab(2)*xij(0);
			RotationRate(1,2) = vab(1)*xij(2)-vab(2)*xij(1);
			RotationRate(1,0) = -RotationRate(0,1);
			RotationRate(2,0) = -RotationRate(0,2);
			RotationRate(2,1) = -RotationRate(1,2);
			RotationRate	  = -0.5 * GK * RotationRate;
			
			// if (StrainRate(2,2)<-1.e-3)
			// cout << "StrainRate"<<StrainRate<<endl;
		
		//}else{ 
			//cout << "Applying kernel corr"<<endl;
			////// New form
			Mat3_t gradv[2],gradvT[2];
			
			//cout<<"gradv"<<gradv[0]<<endl;
			Vec3_t gradK; 
			Mult(GK * P1->gradCorrM,xij,gradK);
			Dyad (vab,gradK,gradv[0]); //outer product. L, velocity gradient tensor
			Mult(GK * P2->gradCorrM,xij,gradK);
			Dyad (vab,gradK,gradv[1]); //outer product. L, velocity gradient tensor
			
			for (int i=0;i<2;i++){
				Trans(gradv[i],gradvT[i]);
				StrainRate_c[i] 	= -0.5*(gradv[i] + gradvT[i]);
				RotationRate_c[i] = -0.5*(gradv[i] - gradvT[i]);
			}
      
      //m_forces_tensors_time += (double)(clock() - m_clock_begin) / CLOCKS_PER_SEC;
			// if (StrainRate_c[0](2,2)<-1.E-3)	
			// cout << "StrainRate_c 1"<<StrainRate_c[0]<<"StrainRate_c 2"<<StrainRate_c[1]<<endl;
			/////////////////////////////////////////////////////////////////////////////////
			///////Chen eq 16. f'xi = Sum_j (mj/dj (fj-fi) Wij,x) / (Sum_j mj/dj (xj-xi) Wij) 
		//}

		// XSPH Monaghan
		if (XSPH != 0.0  && (P1->IsFree*P2->IsFree)) {
			omp_set_lock(&P1->my_lock);
			P1->VXSPH += XSPH*mj/(0.5*(di+dj))*K*-vij;
			omp_unset_lock(&P1->my_lock);

			omp_set_lock(&P2->my_lock);
			P2->VXSPH += XSPH*mi/(0.5*(di+dj))*K*vij;
			omp_unset_lock(&P2->my_lock);
		}

		// Calculating the forces for the particle 1 & 2
		Vec3_t temp = 0.0;
		double temp1 = 0.0;
		Vec3_t temp_c[2];
		double temp1_c[2];
		Vec3_t vc[2];
		
		if (gradKernelCorr){
			for (int i=0;i<2;i++){
				Mult (GKc[i], xij, vc[i]);
			}
		}

		// Original
		// if (GradientType == 0)
			// Mult( GK*xij , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , temp);
		// else
			// Mult( GK*xij , ( 1.0/(di*dj)*(Sigmai + Sigmaj)           + PIij + TIij ) , temp);
    
    clock_begin = clock();
		// NEW
		if (!gradKernelCorr) {
		if (GradientType == 0)
			Mult( GK*xij , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , temp);
		else
			Mult( GK*xij , ( 1.0/(di*dj)*(Sigmai + Sigmaj)           + PIij + TIij ) , temp);
		} else {
				//Should be replaced  dot( xij , GK*xij ) by dot( xij , v )
				//Left in vector form and multiply after??
				for (int i=0;i<2;i++){
					Mult( vc[i] , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , temp_c[i]);
				}
		}//Grad Corr
    
    m_forces_momentum_time += (double)(clock() - clock_begin) / CLOCKS_PER_SEC;
    
		// if (abs(temp(0))>1.e-3){
		// cout << "Strain Rate"<<StrainRate<<endl;
		// cout << "GK*xij"<<GK*xij<<endl;
		// cout << "temp"<<temp<<endl;
		// cout << "Strain Rate c 1"<<StrainRate_c[0]<<", StrainRate c 1"<<StrainRate_c[1]<<endl;
		// cout << "GK*xij corr 1 "<<vc[0] <<"GK*xij corr 2 "<<vc[1]<<endl;
		// cout << "temp_c 1"<<temp_c[0]<<", tempc 2 "<<temp_c[1]<<endl;
		// }
		if (Dimension == 2) temp(2) = 0.0;
		
		if (!gradKernelCorr){
			temp1 = dot( vij , GK*xij );
		} else {
			for (int i=0;i<2;i++){			//TODO: DO THIS ONCE!
				temp1_c[i] = dot( vij , vc[i] );
			}
		}
    
    clock_begin = clock();
		// Locking the particle 1 for updating the properties
		//omp_set_lock(&P1->my_lock);
			if (!gradKernelCorr){
				P1->a					+= mj * temp;
				P1->dDensity	+= mj * (di/dj) * temp1;
			} else{
				P1->a					+= mj * temp_c[0];
				P1->dDensity	+= mj * (di/dj) * temp1_c[0];
			}
			
			

				
			if (P1->IsFree) {
				float mj_dj= mj/dj;
				P1->ZWab	+= mj_dj* K;
				if (!gradKernelCorr){
					P1->StrainRate 		= P1->StrainRate + mj_dj*StrainRate;
					P1->RotationRate 	= P1->RotationRate + mj_dj*RotationRate;
				}
				else {
					P1->StrainRate 		= P1->StrainRate 		+ mj_dj * StrainRate_c[0];
					P1->RotationRate 	= P1->RotationRate 	+ mj_dj * RotationRate_c[0];
				}

			}
			else
				P1->ZWab	= 1.0;

			if (P1->Shepard)
				if (P1->ShepardCounter == P1->ShepardStep)
					P1->SumDen += mj*    K;
		//omp_unset_lock(&P1->my_lock);

		// Locking the particle 2 for updating the properties
		//omp_set_lock(&P2->my_lock);
			if (!gradKernelCorr){
				P2->a					-= mi * temp;
				P2->dDensity	+= mi * (dj/di) * temp1;							
			}else {
				P2->a					-= mi * temp_c[1];
				P2->dDensity	+= mi * (dj/di) * temp1_c[1];
			}
	

			if (P2->IsFree) {
				float mi_di = mi/di;
				P2->ZWab	+= mi_di* K;
				if (!gradKernelCorr){
					P2->StrainRate	 = P2->StrainRate + mi_di*StrainRate;
					P2->RotationRate = P2->RotationRate + mi_di*RotationRate;
				} else {
					P2->StrainRate = P2->StrainRate 		+ mi_di*StrainRate_c[1];
					P2->RotationRate = P2->RotationRate + mi_di*RotationRate_c[1];
				}


			}
			else
				P2->ZWab	= 1.0;

			if (P2->Shepard)
				if (P2->ShepardCounter == P2->ShepardStep)
					P2->SumDen += mi*    K;

		//omp_unset_lock(&P2->my_lock);
 
		//omp_set_lock(&dom_lock); //THIS CAUSES EXTREMELY LONG TIMES
    m_forces_update_time += (double)(clock() - clock_begin) / CLOCKS_PER_SEC;
    //omp_unset_lock(&dom_lock);
  //}//Interaction
} 

  };//SPH