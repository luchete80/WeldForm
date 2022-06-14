#include "Domain.h"

namespace SPH{
inline void Domain::InitImplicitSolver(){
  #pragma omp parallel for num_threads(Nproc)
  #ifdef __GNUC__
  for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
  #else
  for (int i=0; i<Particles.Size(); i++)//Like in Domain::Move
  #endif
  {
    Matrix *mat = new Matrix(6,3);
    Particles[i]->m_B = *mat;
  }
}
inline void Domain::CalcBMat  (){
  Particle *P1, *P2;
  
	#pragma omp parallel for schedule (static) private (P1,P2) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t k=0; k<Nproc;k++) 
	#else
	for (int k=0; k<Nproc;k++) 
	#endif	
	{
    for (size_t p=0; p<SMPairs[k].Size();p++) {
      P1	= Particles[SMPairs[k][p].first];
      P2	= Particles[SMPairs[k][p].second];	 
      double di=0.0,dj=0.0,mi=0.0,mj=0.0;      
			di = P1->Density;
			mi = P1->Mass;

      double h	= (P1->h+P2->h)/2;
      Vec3_t xij	= P1->x - P2->x;
      double rij	= norm(xij);

			dj = P2->Density;
			mj = P2->Mass;    
      double GK	= GradKernel(Dimension, KernelType, rij/h, h);
      Vec3_t d_dx = GK * xij;
      
      omp_set_lock(&P1->my_lock);
      //Like Fraser 11 22 33 23 13 12
      for (int i=0;i<3;i++)
        P1->m_B.Add(i,i, d_dx(i));
      P1->m_B.Add(3,1, d_dx(2));        P1->m_B.Add(3,2, d_dx(1));
      P1->m_B.Add(4,0, d_dx(2));        P1->m_B.Add(4,1, d_dx(2));
      P1->m_B.Add(5,0, d_dx(2));        P1->m_B.Add(5,1, d_dx(0));
      omp_unset_lock(&P1->my_lock);
    }//Pairs
    
  }//For pairs
  
  
}

};