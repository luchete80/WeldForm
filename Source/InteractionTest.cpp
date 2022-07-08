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

};//SPH