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

#include "Particle.h"

namespace SPH {

inline Particle::Particle(int Tag, Vec3_t const & x0, Vec3_t const & v0, double Mass0, double Density0, double h0,bool Fixed)
{
	ct = 0;
	a = 0.0;
    x = x0;
    n = 0.0;
    n0 = 0.0;
    k = 0.0;
    k2 = 0.0;

    Cs		= 0.0;
    P0		= 0.0;
    PresEq	= 0;
    Alpha	= 0.0;
    Beta	= 0.0;

    va = 0.0;
    vb = 0.0;
    NSv = 0.0;
    FSINSv = 0.0;
    v = v0;
    VXSPH = 0.0;
    TI		= 0.0;
    TIn		= 4.0;
    TIInitDist  = 0.0;

    Densitya = 0.0;
    Densityb = 0.0;
    Density = Density0;
    RefDensity = Density0;

    Mass = Mass0;
    FPMassC = 1.0;
    IsFree = !Fixed;
    h = h0;
    Pressure=0.0;
    FSIPressure=0.0;
    ID = Tag;
    CC[0]= CC[1] = CC[2] = 0;
    LL=0;
    ZWab = 0.0;
    SumDen = 0.0;
    dDensity=0.0;
    ShearRate = 0.0;
    MuRef = Mu = 0.0;
		VisM = 0;
    T0 = 0.0;
    m = 300.0;
    SumKernel = 0.0;
    FSISumKernel = 0.0;
    G = 0.0;
    K = 0.0;
    Material = 0;
    Fail = 0;
    c = 0.0;
    phi = 0.0;
    psi = 0.0;
    d =0.0;
    Sigmay = 0.0;
    NoSlip = false;
    Shepard = false;
    InOut = 0;
    FirstStep = true;
    V = Mass/RefDensity;
    RhoF = 0.0;
    IsSat = false;
    SatCheck = false;
    ShepardStep = 40;
    ShepardCounter = 0;
    S = 0.0;
    VarPorosity = false;
    SeepageType = 0;
    S = 0;
	LES = false;
	SBar = 0.0;
	CSmag = 0.17;



    set_to_zero(Strainb);
    set_to_zero(Strain);
    set_to_zero(Sigmab);
    set_to_zero(Sigma);
//    set_to_zero(FSISigma);
    set_to_zero(Sigmaa);
    set_to_zero(ShearStress);
    set_to_zero(ShearStressb);
    set_to_zero(TIR);
    set_to_zero(StrainRate);
    set_to_zero(RotationRate);
    omp_init_lock(&my_lock);

}

inline void Particle::Move(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin, size_t Scheme, Mat3_t I)
{
	if (Scheme == 0)
		Move_MVerlet(I, dt);
	else
		Move_Leapfrog(I, dt);


	//Periodic BC particle position update
	if (Domainsize(0)>0.0)
	{
		(x(0)>(domainmax(0))) ? x(0) -= Domainsize(0) : x(0);
		(x(0)<(domainmin(0))) ? x(0) += Domainsize(0) : x(0);
	}
	if (Domainsize(1)>0.0)
	{
		(x(1)>(domainmax(1))) ? x(1) -= Domainsize(1) : x(1);
		(x(1)<(domainmin(1))) ? x(1) += Domainsize(1) : x(1);
	}
	if (Domainsize(2)>0.0)
	{
		(x(2)>(domainmax(2))) ? x(2) -= Domainsize(2) : x(2);
		(x(2)<(domainmin(2))) ? x(2) += Domainsize(2) : x(2);
	}

}

inline void Particle::Mat1(double dt)
{
	Pressure 	= EOS(PresEq, Cs, P0,Density, RefDensity);
	double temp	= (StrainRate(0,0)*StrainRate(0,0) + 2.0*StrainRate(0,1)*StrainRate(1,0) +
			2.0*StrainRate(0,2)*StrainRate(2,0) + StrainRate(1,1)*StrainRate(1,1) +
			2.0*StrainRate(1,2)*StrainRate(2,1) + StrainRate(2,2)*StrainRate(2,2));

	ShearRate	= sqrt(0.5*temp);
	SBar		= sqrt(2.0*temp);

	// LES model
	if (LES)
	{
		Mu	= MuRef + RefDensity*pow((CSmag*h),2.0)*SBar;
	}

	// Bingham viscosity calculation
	if (T0>0.0)
	{
		switch (VisM)
		{
			case 0:
			// Bingham
				if (ShearRate !=0.0)
					Mu = MuRef + T0*(1-exp(-m*ShearRate))/ShearRate;
				else
					Mu = MuRef + T0*m;
				break;
			case 1:
			// Cross
				Mu = (1000.0*MuRef + MuRef*MuRef*1000.0/T0*ShearRate)/(1+1000.0*MuRef/T0*ShearRate);
				break;
			default:
				std::cout << "Non-Newtonian Viscosity Type No is out of range. Please correct it and run again" << std::endl;
				std::cout << "0 => Bingham" << std::endl;
				std::cout << "1 => Cross" << std::endl;
				abort();
				break;
		}
	}
}


inline void Particle::Move_MVerlet (Mat3_t I, double dt)
{
	if (FirstStep)
	{
		ct = 30;
		FirstStep = false;
	}

	x += dt*(v+VXSPH) + 0.5*dt*dt*a;

	if (ct == 30)
	{
		if (Shepard && ShepardCounter == ShepardStep)
		{
			if (ZWab>0.6)
			{
				Densityb	= SumDen/ZWab;
//				Densityb	= Density;
				Density		= SumDen/ZWab;
			}
			else
			{
				Densityb	= Density;
				Density		+=dt*dDensity;
			}
		}
		else
		{
			Densityb		= Density;
			Density			+=dt*dDensity;
		}

		vb	= v;
		v	+=dt*a;
	}
	else
	{
		if (Shepard && ShepardCounter == ShepardStep)
		{
			if (ZWab>0.6)
			{
				Densityb	= SumDen/ZWab;
//				Densityb	= Density;
				Density		= SumDen/ZWab;
			}
			else
			{
				double dens	= Density;
				Density		= Densityb + 2.0*dt*dDensity;
				Densityb	= dens;
			}
		}
		else
		{
			double dens	= Density;
			Density		= Densityb + 2.0*dt*dDensity;
			Densityb	= dens;
		}

		Vec3_t temp;
		temp	= v;
		v		= vb + 2*dt*a;
		vb		= temp;
	}

	switch (Material)
    {case 1:
    	Mat1(dt);
		break;
    case 2:
    	Mat2MVerlet(dt);
    	break;
    // case 3:
    	// Mat3MVerlet(I,dt);
    	// break;
   default:
	   	std::cout << "Material Type No is out of range. Please correct it and run again" << std::endl;
		std::cout << "1 => Fluid" << std::endl;
		std::cout << "2 => Solid" << std::endl;
		std::cout << "3 => Soil" << std::endl;
	    abort();
	    break;
    }
	if (ct == 30) ct = 0; else ct++;
	if (ShepardCounter == ShepardStep) ShepardCounter = 0; else ShepardCounter++;
}

inline void Particle::Mat2MVerlet(double dt)
{
	Pressure = EOS(PresEq, Cs, P0,Density, RefDensity);

	// Jaumann rate terms
	Mat3_t RotationRateT, Stress,SRT,RS;
	Trans(RotationRate,RotationRateT);
	Mult(ShearStress,RotationRateT,SRT);
	Mult(RotationRate,ShearStress,RS);

	// Elastic prediction step (ShearStress_e n+1)
	Stress			= ShearStress;
	if (ct == 30)
		ShearStress	= dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStress;
	else
		ShearStress	= 2.0*dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStressb;
	ShearStressb	= Stress;

	if (Fail == 1)
	{
		double J2	= 0.5*(ShearStress(0,0)*ShearStress(0,0) + 2.0*ShearStress(0,1)*ShearStress(1,0) +
						2.0*ShearStress(0,2)*ShearStress(2,0) + ShearStress(1,1)*ShearStress(1,1) +
						2.0*ShearStress(1,2)*ShearStress(2,1) + ShearStress(2,2)*ShearStress(2,2));
		//Scale back
		ShearStress	= std::min((Sigmay/sqrt(3.0*J2)),1.0)*ShearStress;
	}

	Sigma			= -Pressure * OrthoSys::I + ShearStress;

	Stress	= Strain;
	if (ct == 30)
		Strain	= dt*StrainRate + Strain;
	else
		Strain	= 2.0*dt*StrainRate + Strainb;
	Strainb	= Stress;


	if (Fail > 1)
	{
		std::cout<<"Undefined failure criteria for solids"<<std::endl;
		abort();
	}
}

inline void Particle::Move_Leapfrog(Mat3_t I, double dt)
{
	if (FirstStep)
	{
		Densitya = Density - dt/2.0*dDensity;
		va = v - dt/2.0*a;
	}
	Densityb = Densitya;
	Densitya += dt*dDensity;
	Density = (Densitya+Densityb)/2.0;
	vb = va;
	va += dt*a;
	v = (va + vb)/2.0;
	x += dt*va;

	switch (Material)
    {case 1:
    	Mat1(dt);
		break;
    case 2:
    	Mat2Leapfrog(dt);
    	break;
    // case 3:
    	// Mat3Leapfrog(I,dt);
    	// break;
   default:
	   	std::cout << "Material Type No is out of range. Please correct it and run again" << std::endl;
		std::cout << "1 => Fluid" << std::endl;
		std::cout << "2 => Solid" << std::endl;
		std::cout << "3 => Soil" << std::endl;
	    abort();
	    break;
    }
	if (FirstStep) FirstStep = false;

}

inline void Particle::Mat2Leapfrog(double dt)
{
	Pressure = EOS(PresEq, Cs, P0,Density, RefDensity);

	// Jaumann rate terms
	Mat3_t RotationRateT,SRT,RS;
	Trans(RotationRate,RotationRateT);
	Mult(ShearStress,RotationRateT,SRT);
	Mult(RotationRate,ShearStress,RS);

	// Elastic prediction step (ShearStress_e n+1)
	if (FirstStep)
		ShearStressa	= -dt/2.0*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStress;

	ShearStressb	= ShearStressa;
	ShearStressa	= dt*(2.0*G*(StrainRate-1.0/3.0*(StrainRate(0,0)+StrainRate(1,1)+StrainRate(2,2))*OrthoSys::I)+SRT+RS) + ShearStressa;

	if (Fail == 1)
	{
		double J2	= 0.5*(ShearStressa(0,0)*ShearStressa(0,0) + 2.0*ShearStressa(0,1)*ShearStressa(1,0) +
						2.0*ShearStressa(0,2)*ShearStressa(2,0) + ShearStressa(1,1)*ShearStressa(1,1) +
						2.0*ShearStressa(1,2)*ShearStressa(2,1) + ShearStressa(2,2)*ShearStressa(2,2));
		//Scale back
		ShearStressa= std::min((Sigmay/sqrt(3.0*J2)),1.0)*ShearStressa;
	}
	ShearStress	= 1.0/2.0*(ShearStressa+ShearStressb);

	Sigma = -Pressure * OrthoSys::I + ShearStress;

	if (FirstStep)
		Straina	= -dt/2.0*StrainRate + Strain;
	Strainb	= Straina;
	Straina	= dt*StrainRate + Straina;
	Strain	= 1.0/2.0*(Straina+Strainb);


	if (Fail > 1)
	{
		std::cout<<"Undefined failure criteria for solids"<<std::endl;
		abort();
	}
}



inline void Particle::translate(double dt, Vec3_t Domainsize, Vec3_t domainmax, Vec3_t domainmin)
{
	x = x + dt*v + 0.5*dt*dt*a;

	// Evolve velocity
	Vec3_t temp;
	temp = v;
	v = vb + 2*dt*a;
	vb = temp;

	//Periodic BC particle position update
	if (Domainsize(0)>0.0)
	{
		(x(0)>(domainmax(0))) ? x(0) -= Domainsize(0) : x(0);
		(x(0)<(domainmin(0))) ? x(0) += Domainsize(0) : x(0);
	}
	if (Domainsize(1)>0.0)
	{
		(x(1)>(domainmax(1))) ? x(1) -= Domainsize(1) : x(1);
		(x(1)<(domainmin(1))) ? x(1) += Domainsize(1) : x(1);
	}
	if (Domainsize(2)>0.0)
	{
		(x(2)>(domainmax(2))) ? x(2) -= Domainsize(2) : x(2);
		(x(2)<(domainmin(2))) ? x(2) += Domainsize(2) : x(2);
	}
}

}; // namespace SPH
