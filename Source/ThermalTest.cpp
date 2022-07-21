#include "Domain.h"
#include <vector>
using namespace std;

namespace SPH {


inline double Domain::CalcTempIncPair (Particle *P1, Particle *P2) {
	double di=0.0,dj=0.0,mi=0.0,mj=0.0;
	
	std::vector < double> temp(Particles.Size());
	
  double ret =0.;

  Vec3_t xij;

  xij	= P1->x - P2->x;
  double h	= (P1->h+P2->h)/2.0;
  double GK	= GradKernel(Dimension, KernelType, norm(xij)/h, h);	    
  
  di = P1->Density; mi = P1->Mass;
  dj = P2->Density; mj = P2->Mass;
  
  //ret = mj/dj * 4. * ( P1->k_T * P2->k_T) / (P1->k_T + P2->k_T) * ( P1->T - P2->T) * dot( xij , GK*xij )/ (norm(xij)*norm(xij));

	return ret;
}

//THIS IS REDUNDAND, LIKE FRASER ALGORITHM ,CALCULATING EACH PAIR 2 TIMES
inline void Domain::CalcTempIncPP() {
  Particle *P1, *P2;
  double k;
	#pragma omp parallel for schedule (static) private (P1,P2,k) num_threads(Nproc)
  for (int i=0; i < Particles.Size();i++){
    P1	= Particles[i]; 
    P1->dTdt = 0.;
    k = 1./(P1->Density * P1->cp_T);
    for (int n=0;n<ipair_SM[i];n++){      
      P2 = Particles[Anei[i][n]];
      P1->dTdt += CalcTempIncPair(P1,P2);
    }
    for (int n=0;n<jpair_SM[i];n++) {
      P2 = Particles[Anei[i][MAX_NB_PER_PART-1-n]];
      P1->dTdt += CalcTempIncPair(P1,P2);
    }
    P1->dTdt *= k; 
  }
}

};//sph