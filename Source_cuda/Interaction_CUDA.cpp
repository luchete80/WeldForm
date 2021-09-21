/***********************************************************************************
* PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids        *
*             and soils) using Smoothed Particle Hydrodynamics method              *
* Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo-Torres              *
*                                                                                  *
* This file is part of PersianSPH                                                  *
*                                                                                  *
* This is free software; you can redistribute it and/or modify it under the        *
* terms of the GNU General Public License as published by the Free Software        *
* Foundation; either version 3 of the License, or (at your option) any later       *
* version.                                                                         *
*                                                                                  *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY  *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.         *
*                                                                                  *
* You should have received a copy of the GNU General Public License along with     *
* PersianSPH; if not, see <http://www.gnu.org/licenses/>                           *
************************************************************************************/

#include "Domain.h"

namespace SPH {

#ifdef CUDA

class Particle_cu
{
    public :
        int value;
        float rate;
		Vec3_t	x;		///< Position of the particle n
		Vec3_t	vb;		///< Velocity of the particle n-1 (Modified Verlet)
		Vec3_t	va;		///< Velocity of the particle n+1/2 (Leapfrog)
		Vec3_t	v;		///< Velocity of the particle n+1
		Vec3_t	NSv;		///< Velocity of the fixed particle for no-slip BC
		Vec3_t	VXSPH;		///< Mean Velocity of neighbor particles for updating the particle position (XSPH)
		Vec3_t	a;		///< Acceleration of the particle n
		
        __device__ __host__ Particle_cu()
        {
            value = 0; rate = 0;
        }
        __device__ __host__ Particle_cu(int v,float r)
        {
            value = v; rate = r;
        }
        __device__ __host__ ~Particle_cu() {};
}

__global__
inline void Domain::CalcForce2233_CUDA(Particle_cu * P1, Particle_cu * P2)
{
	double h	= (P1->h+P2->h)/2;
	Vec3_t xij	= P1->x - P2->x;

	Periodic_X_Correction(xij, h, P1, P2);

	double rij	= norm(xij);

	if ((rij/h)<=Cellfac)
	{
		double di=0.0,dj=0.0,mi=0.0,mj=0.0;
		double Alpha	= (P1->Alpha + P2->Alpha)/2.0;
		double Beta	= (P1->Beta + P2->Beta)/2.0;


		if (!P1->IsFree) {
			di = DensitySolid(P2->PresEq, P2->Cs, P2->P0,P1->Pressure, P2->RefDensity);
			mi = P1->FPMassC * P2->Mass;
		} else {
			di = P1->Density;
			mi = P1->Mass;
		}
		if (!P2->IsFree) {
			dj = DensitySolid(P1->PresEq, P1->Cs, P1->P0,P2->Pressure, P1->RefDensity);
			mj = P2->FPMassC * P1->Mass;
		} else {
			dj = P2->Density;
			mj = P2->Mass;
		}
		

		Vec3_t vij	= P1->v - P2->v;
		
		double GK	= GradKernel(Dimension, KernelType, rij/h, h);
		double K	= Kernel(Dimension, KernelType, rij/h, h);
		
		// double GK	= m_kernel.gradW(rij/h);
		// double K	= m_kernel.W(rij/h);
		
		// Artificial Viscosity
		Mat3_t PIij;
		set_to_zero(PIij);
		if (Alpha!=0.0 || Beta!=0.0)
		{
			double MUij = h*dot(vij,xij)/(rij*rij+0.01*h*h);					///<(2.75) Li, Liu Book
			double Cij;
			double Ci,Cj;
			if (!P1->IsFree) Ci = SoundSpeed(P2->PresEq, P2->Cs, di, P2->RefDensity); else Ci = SoundSpeed(P1->PresEq, P1->Cs, di, P1->RefDensity);
			if (!P2->IsFree) Cj = SoundSpeed(P1->PresEq, P1->Cs, dj, P1->RefDensity); else Cj = SoundSpeed(P2->PresEq, P2->Cs, dj, P2->RefDensity);
			Cij = 0.5*(Ci+Cj);
			
			if (dot(vij,xij)<0) PIij = (Alpha*Cij*MUij+Beta*MUij*MUij)/(0.5*(di+dj)) * I;		///<(2.74) Li, Liu Book
		}

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

		// Calculation strain rate tensor
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

		// Calculation rotation rate tensor
		RotationRate(0,1) = vab(0)*xij(1)-vab(1)*xij(0);
		RotationRate(0,2) = vab(0)*xij(2)-vab(2)*xij(0);
		RotationRate(1,2) = vab(1)*xij(2)-vab(2)*xij(1);
		RotationRate(1,0) = -RotationRate(0,1);
		RotationRate(2,0) = -RotationRate(0,2);
		RotationRate(2,1) = -RotationRate(1,2);
		RotationRate	  = -0.5 * GK * RotationRate;

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

		if (GradientType == 0)
			Mult( GK*xij , ( 1.0/(di*di)*Sigmai + 1.0/(dj*dj)*Sigmaj + PIij + TIij ) , temp);
		else
			Mult( GK*xij , ( 1.0/(di*dj)*(Sigmai + Sigmaj)           + PIij + TIij ) , temp);

		if (Dimension == 2) temp(2) = 0.0;
		temp1 = dot( vij , GK*xij );

		// Locking the particle 1 for updating the properties
		//omp_set_lock(&P1->my_lock);
			P1->a		+= mj * temp;
			P1->dDensity	+= mj * (di/dj) * temp1;

			if (P1->IsFree) {
				float mj_dj= mj/dj;
				P1->ZWab	+= mj_dj* K;
				P1->StrainRate	 = P1->StrainRate + mj_dj*StrainRate;
				P1->RotationRate = P1->RotationRate + mj_dj*RotationRate;
			}
			else
				P1->ZWab	= 1.0;

			if (P1->Shepard)
				if (P1->ShepardCounter == P1->ShepardStep)
					P1->SumDen += mj*    K;
		//omp_unset_lock(&P1->my_lock);

		// Locking the particle 2 for updating the properties
		//omp_set_lock(&P2->my_lock);
			P2->a		-= mi * temp;
			P2->dDensity	+= mi * (dj/di) * temp1;
			if (P2->IsFree) {
				float mi_di = mi/di;
				P2->ZWab	+= mi_di* K;
				P2->StrainRate	 = P2->StrainRate + mi_di*StrainRate;
				P2->RotationRate = P2->RotationRate + mi_di*RotationRate;

			}
			else
				P2->ZWab	= 1.0;

			if (P2->Shepard)
				if (P2->ShepardCounter == P2->ShepardStep)
					P2->SumDen += mi*    K;

		//omp_unset_lock(&P2->my_lock);
	}
}

#endif


}; // namespace SPH
