#include "Domain.h"
#include <chrono>
//#include <time.h>       /* time_t, struct tm, difftime, time, mktime */
#include <ctime> //Clock
#include <vector>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <iostream>
#include <numeric>


using namespace std;

namespace SPH {

__global__	
inline void Domain::PrimaryComputeAcceleration_CUDA () {
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t k=0; k<Nproc;k++) 
	#else
	for (int k=0; k<Nproc;k++) 
	#endif
	
	{
		size_t P1,P2;
		Vec3_t xij;
		double h,K;
		// Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
		for (size_t a=0; a<FSMPairs[k].Size();a++) {
			P1	= FSMPairs[k][a].first;
			P2	= FSMPairs[k][a].second;
			xij	= Particles[P1]->x-Particles[P2]->x;
			h	= (Particles[P1]->h+Particles[P2]->h)/2.0;

			Periodic_X_Correction(xij, h, Particles[P1], Particles[P2]);

			K	= Kernel(Dimension, KernelType, norm(xij)/h, h);
			if ( !Particles[P1]->IsFree ) {
				omp_set_lock(&Particles[P1]->my_lock);
										Particles[P1]->SumKernel+= K;
					Particles[P1]->Pressure	+= Particles[P2]->Pressure * K + dot(Gravity,xij)*Particles[P2]->Density*K;
					Particles[P1]->Sigma 	 = Particles[P1]->Sigma + K * Particles[P2]->Sigma;
					if (Particles[P1]->NoSlip)		Particles[P1]->NSv 	+= Particles[P2]->v * K;
				omp_unset_lock(&Particles[P1]->my_lock);
			} else {
				omp_set_lock(&Particles[P2]->my_lock);
										Particles[P2]->SumKernel+= K;
					Particles[P2]->Pressure	+= Particles[P1]->Pressure * K + dot(Gravity,xij)*Particles[P1]->Density*K;
					Particles[P2]->Sigma	 = Particles[P2]->Sigma + K * Particles[P1]->Sigma;
					if (Particles[P2]->NoSlip)		Particles[P2]->NSv 	+= Particles[P1]->v * K;
				omp_unset_lock(&Particles[P2]->my_lock);
			}
		}
	}

	// Calculateing the finala value of the smoothed pressure, velocity and stress for fixed particles
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t i=0; i<FixedParticles.Size(); i++)
	#else
	for (int i=0; i<FixedParticles.Size(); i++)
	#endif

		if (Particles[FixedParticles[i]]->SumKernel!= 0.0) {
			size_t a = FixedParticles[i];
			Particles[a]->Pressure	= Particles[a]->Pressure/Particles[a]->SumKernel;
			Particles[a]->Sigma	= 1.0/Particles[a]->SumKernel*Particles[a]->Sigma;
			if (Particles[a]->NoSlip)	Particles[a]->NSv	= Particles[a]->NSv/Particles[a]->SumKernel;

			// Tensile Instability for fixed soil and solid particles
			if (Particles[a]->TI > 0.0)
			{
				// XY plane must be used, It is very slow in 3D
				if (Dimension == 2) {
					double teta, Sigmaxx, Sigmayy, C, S;
					if ((Particles[a]->Sigma(0,0)-Particles[a]->Sigma(1,1))!=0.0)
						teta = 0.5*atan(2.0*Particles[a]->Sigma(0,1)/(Particles[a]->Sigma(0,0)-Particles[a]->Sigma(1,1)));
					else
						teta = M_PI/4.0;

					C = cos(teta);
					S = sin(teta);
					Sigmaxx = C*C*Particles[a]->Sigma(0,0) + 2.0*C*S*Particles[a]->Sigma(0,1) + S*S*Particles[a]->Sigma(1,1);
					Sigmayy = S*S*Particles[a]->Sigma(0,0) - 2.0*C*S*Particles[a]->Sigma(0,1) + C*C*Particles[a]->Sigma(1,1);
					if (Sigmaxx>0) Sigmaxx = -Particles[a]->TI * Sigmaxx/(Particles[a]->Density*Particles[a]->Density); else Sigmaxx = 0.0;
					if (Sigmayy>0) Sigmayy = -Particles[a]->TI * Sigmayy/(Particles[a]->Density*Particles[a]->Density); else Sigmayy = 0.0;
					Particles[a]->TIR(0,0) = C*C*Sigmaxx + S*S*Sigmayy;
					Particles[a]->TIR(1,1) = S*S*Sigmaxx + C*C*Sigmayy;
					Particles[a]->TIR(0,1) = Particles[a]->TIR(1,0) = S*C*(Sigmaxx-Sigmayy);
				}
				else {
					Mat3_t Vec,Val,VecT,temp;
					Rotation(Particles[a]->Sigma,Vec,VecT,Val);
					double pc_ti_inv_d2=Particles[a]->TI/(Particles[a]->Density*Particles[a]->Density);//Precompute some values
					// if (Val(0,0)>0) Val(0,0) = -Particles[a]->TI * Val(0,0)/(Particles[a]->Density*Particles[a]->Density); else Val(0,0) = 0.0;
					// if (Val(1,1)>0) Val(1,1) = -Particles[a]->TI * Val(1,1)/(Particles[a]->Density*Particles[a]->Density); else Val(1,1) = 0.0;
					// if (Val(2,2)>0) Val(2,2) = -Particles[a]->TI * Val(2,2)/(Particles[a]->Density*Particles[a]->Density); else Val(2,2) = 0.0;
					if (Val(0,0)>0) Val(0,0) = -pc_ti_inv_d2 * Val(0,0); else Val(0,0) = 0.0;
					if (Val(1,1)>0) Val(1,1) = -pc_ti_inv_d2 * Val(1,1); else Val(1,1) = 0.0;
					if (Val(2,2)>0) Val(2,2) = -pc_ti_inv_d2 * Val(2,2); else Val(2,2) = 0.0;

					Mult(Vec,Val,temp);
					Mult(temp,VecT,Particles[a]->TIR);
				}
			}
		}


}

}; // namespace SPH