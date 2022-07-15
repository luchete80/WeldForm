#include "Domain.h"

///////////////////////////////////////////////////////////////////
// ONLY FOR DELETION: IS FOR TESTING PARALLELIZE BY PARTICLE OR BY PAIR
// PARALLELIZATION PER PARTICLE
#define ID_TEST 0
namespace SPH{
// ONLY FOR TESTING  
inline void Domain::CalcAccelPair(Particle * P1, Particle * P2){
	double h	= (P1->h+P2->h)/2;
	Vec3_t xij	= P1->x - P2->x;

	Periodic_X_Correction(xij, h, P1, P2);

	double rij	= norm(xij);

  double clock_begin;
  Vec3_t temp = 0.;
	// if ((rij/h)<=Cellfac)
	// {
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
		double Alpha	= (P1->Alpha + P2->Alpha)/2.0;
		double Beta	= (P1->Beta + P2->Beta)/2.0;

    di = P1->Density;
    mi = P1->Mass;
    dj = P2->Density;
    mj = P2->Mass;

		Vec3_t vij	= P1->v - P2->v;
		
		double GK	= GradKernel(Dimension, KernelType, rij/h, h);
		double K	= Kernel(Dimension, KernelType, rij/h, h);
		

		// Artificial Viscosity
		Mat3_t PIij;
		set_to_zero(PIij);
		if (Alpha!=0.0 || Beta!=0.0) {
			double MUij = h*dot(vij,xij)/(rij*rij+0.01*h*h);					///<(2.75) Li, Liu Book
			double Cij;
			double Ci,Cj;
			if (!P1->IsFree) Ci = SoundSpeed(P2->PresEq, P2->Cs, di, P2->RefDensity); else Ci = SoundSpeed(P1->PresEq, P1->Cs, di, P1->RefDensity);
			if (!P2->IsFree) Cj = SoundSpeed(P1->PresEq, P1->Cs, dj, P1->RefDensity); else Cj = SoundSpeed(P2->PresEq, P2->Cs, dj, P2->RefDensity);
			Cij = 0.5*(Ci+Cj);
			
			if (dot(vij,xij)<0) PIij = (Alpha*Cij*MUij+Beta*MUij*MUij)/(0.5*(di+dj)) * I;		///<(2.74) Li, Liu Book
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
		if (P1->TI > 0.0 || P2->TI > 0.0) 
			TIij = pow((K/Kernel(Dimension, KernelType, (P1->TIInitDist + P2->TIInitDist)/(2.0*h), h)),(P1->TIn+P2->TIn)/2.0)*(P1->TIR+P2->TIR);
			//TIij = pow((K/m_kernel.W((P1->TIInitDist + P2->TIInitDist)/(2.0*h))),(P1->TIn+P2->TIn)/2.0)*(P1->TIR+P2->TIR);

		if (!gradKernelCorr) {
		if (GradientType == 0)
			Mult( GK*xij , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , temp);
    }
    if (!gradKernelCorr){
      //P1->a					+= mj * temp;
			P2->a					-= mi * temp;
    }
			
}

//////////////////////////
//// ONLY FOR TESTING FUNCTION ///
//// PARALLELIZATION BY PARTICLE
inline void Domain::CalcAccelPP() {
  Particle *P1, *P2;

	#pragma omp parallel for schedule (static) private (P1,P2) num_threads(Nproc)
  for (int i=0; i < Particles.Size();i++){
    for (int n=0;n<ipair_SM[i];n++){
      P1	= Particles[i]; 
      P2 = Particles[Anei[i][n]];
      CalcAccelPair(P1,P2);
    }
    for (int n=0;n<jpair_SM[i];n++) {
      P1	= Particles[i]; 
      P2 = Particles[Anei[i][MAX_NB_PER_PART-1-n]];
      CalcAccelPair(P1,P2);
    }
  }
}

// ONLY WRITES IN P2

//Similar but not densities
inline void Domain::CalcRateTensorsPair (Particle *P1, Particle *P2) {          

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
		
		if (!gradKernelCorr)
			temp1 = dot( vij , GK*xij );
    
    clock_begin = clock();
		// Locking the particle 1 for updating the properties

    float mj_dj= mj/dj;

    float mi_di = mi/di;
    P2->StrainRate	 = P2->StrainRate + mi_di*StrainRate;
    P2->RotationRate = P2->RotationRate + mi_di*RotationRate;


}

inline void Domain::CalcTensorsPP() {
  Particle *P1, *P2;

	#pragma omp parallel for schedule (static) private (P1,P2) num_threads(Nproc)
  for (int i=0; i < Particles.Size();i++){
    for (int n=0;n<ipair_SM[i];n++){
      P1	= Particles[i]; 
      P2 = Particles[Anei[i][n]];
      CalcRateTensorsPair(P1,P2);
    }
    for (int n=0;n<jpair_SM[i];n++) {
      P1	= Particles[i]; 
      P2 = Particles[Anei[i][MAX_NB_PER_PART-1-n]];
      CalcRateTensorsPair(P1,P2);
    }
  }
}

inline void Domain::CalcDensIncPairs(Particle *P1, Particle *P2) {

      double h	= (P1->h+P2->h)/2;
      Vec3_t xij	= P1->x - P2->x;
      double rij	= norm(xij);

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
      if (!gradKernelCorr){
        P2->dDensity	+= mi * (dj/di) * temp1;							
      }else {
        P2->dDensity	+= mi * (dj/di) * temp1_c[1];
      }
}

inline void Domain::CalcDensPP() {
  Particle *P1, *P2;

	#pragma omp parallel for schedule (static) private (P1,P2) num_threads(Nproc)
  for (int i=0; i < Particles.Size();i++){
    for (int n=0;n<ipair_SM[i];n++){
      P1	= Particles[i]; 
      P2 = Particles[Anei[i][n]];
      CalcDensIncPairs(P1,P2);
    }
    for (int n=0;n<jpair_SM[i];n++) {
      P1	= Particles[i]; 
      P2 = Particles[Anei[i][MAX_NB_PER_PART-1-n]];
      CalcDensIncPairs(P1,P2);
    }
  }
}

};//SPH