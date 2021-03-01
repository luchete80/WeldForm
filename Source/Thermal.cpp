#include "Domain.h"

namespace SPH {

inline void Domain::CalcTempInc () {
	double di=0.0,dj=0.0,mi=0.0,mj=0.0;
	//#pragma omp parallel for schedule (static) num_threads(Nproc)
	for ( size_t k = 0; k < Nproc ; k++) {
		Particle *P1,*P2;
		Vec3_t xij;
		double h,GK;
		//TODO: DO THE LOCK PARALLEL THING
		// Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
		double temp;
		for (size_t a=0; a<SMPairs[k].Size();a++) {
			P1	= Particles[SMPairs[k][a].first];
			P2	= Particles[SMPairs[k][a].second];
			xij	= Particles[P1]->x-Particles[P2]->x;
			h	= (Particles[P1]->h+Particles[P2]->h)/2.0;
			GK	= GradKernel(Dimension, KernelType, norm(xij)/h, h);
			di = P1->Density; mi = P1->Mass;
			dj = P2->Density; mj = P2->Mass;
			//Frasier  Eqn 3.99 dTi/dt= 1/(rhoi_CPi) * Sum_j(mj/rho_j * 4*ki kj/ (ki + kj ) (Ti - Tj)  ) 
			temp += mj/dj * 4. * ( P1->k_T * P2->k_T) / (P1->k_T + P2->k_T) * ( P1->T * P2->T) * dot( xij , GK*xij );
			omp_set_lock(&P1->my_lock);
				P1->T		+= 1./(di*cp_T) * temp;
			omp_unset_lock(&P1->my_lock);
			// Locking the particle 2 for updating the properties
			omp_set_lock(&P2->my_lock);
				P2->T		-= 1./(di*cp_T) * temp;
			omp_unset_lock(&P2->my_lock);
		}
		// for (size_t i=0; i<FSMPairs[k].Size();i++)
			// CalcForce2233(Particles[FSMPairs[k][i].first],Particles[FSMPairs[k][i].second]);
	}//Nproc

}
inline void Domain::CalcConvHeat(){
	
	//Fraser Eq 3-121
	for ( size_t k = 0; k < Nproc ; k++) {
	
	}		
}
};