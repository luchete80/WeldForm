#include "Domain.h"

namespace SPH {

inline void Domain::InFlowBCFresh()
{
	int temp, temp1;
	int q1,q2,q3;
	if (BC.inoutcounter == 0)
	{
		if (BC.InOutFlow==1 || BC.InOutFlow==3)
		{
			if (!(norm(BC.inv)>0.0) || !(BC.inDensity>0.0))
			{
				std::cout<< "BC.inv or BC.inDensity are not defined, please define them and run the code again."<<std::endl;
				abort();
			}

			BC.InPart.Clear();
			BC.InFlowLoc1  = BLPF(0) + BC.cellfac*hmax;
			temp1 = (int) (floor((BC.InFlowLoc1 - BLPF(0)) / CellSize(0)));

			for (q2=0; BC.Periodic[1]? (q2<(CellNo[1]-2)) : (q2<CellNo[1]) ; q2++)
			for (q3=0; BC.Periodic[2]? (q3<(CellNo[2]-2)) : (q3<CellNo[2]) ; q3++)
			for (q1=0; q1<(temp1 + 1)                                      ; q1++)
			{
				if (HOC[q1][q2][q3]!=-1)
				{
					temp = HOC[q1][q2][q3];
					while (temp != -1)
					{
						if (Particles[temp]->IsFree && (Particles[temp]->x(0) <= BC.InFlowLoc1) )
						{
							BC.InPart.Push(temp);
							Particles[temp]->InOut = 1;
						}
						temp = Particles[temp]->LL;
					}
				}
			}
			BC.InFlowLoc2  = Particles[BC.InPart[0]]->x(0);
			BC.InFlowLoc3  = Particles[BC.InPart[0]]->x(0);
			#pragma omp parallel for schedule(static) num_threads(Nproc)
			for (int i=0 ; i<BC.InPart.Size() ; i++)
			{
				if (Particles[BC.InPart[i]]->x(0) < BC.InFlowLoc2) BC.InFlowLoc2  = Particles[BC.InPart[i]]->x(0);
				if (Particles[BC.InPart[i]]->x(0) > BC.InFlowLoc3) BC.InFlowLoc3  = Particles[BC.InPart[i]]->x(0);
			}
		}

		if (BC.InOutFlow==2 || BC.InOutFlow==3)
			BC.OutFlowLoc = TRPR(0) - BC.cellfac*hmax;

		BC.inoutcounter = 2;
	}

	if (BC.inoutcounter == 1)
	{
		BC.InPart.Clear();
		temp1 = (int) (floor((BC.InFlowLoc1 - BLPF(0)) / CellSize(0)));

		for (q2=0; BC.Periodic[1]? (q2<(CellNo[1]-2)) : (q2<CellNo[1]) ; q2++)
		for (q3=0; BC.Periodic[2]? (q3<(CellNo[2]-2)) : (q3<CellNo[2]) ; q3++)
		for (q1=0; q1<(temp1 + 1)                                      ; q1++)
		{
			if (HOC[q1][q2][q3]!=-1)
			{
				temp = HOC[q1][q2][q3];
				while (temp != -1)
				{
					if (Particles[temp]->IsFree && (Particles[temp]->x(0) <= BC.InFlowLoc1) && Particles[temp]->InOut==1)
						BC.InPart.Push(temp);
					temp = Particles[temp]->LL;
				}
			}
		}
		BC.inoutcounter = 2;
	}


	if (BC.InOutFlow==2 || BC.InOutFlow==3)
	{
		BC.OutPart.Clear();
		temp1 = (int) (floor((BC.OutFlowLoc - BLPF(0)) / CellSize(0)));

		for (q2=0     ; BC.Periodic[1]? (q2<(CellNo[1]-2)) : (q2<CellNo[1]) ; q2++)
		for (q3=0     ; BC.Periodic[2]? (q3<(CellNo[2]-2)) : (q3<CellNo[2]) ; q3++)
		for (q1=temp1 ; q1<CellNo[0]                                        ; q1++)
		{
			if (HOC[q1][q2][q3]!=-1)
			{
				temp = HOC[q1][q2][q3];
				while (temp != -1)
				{
					if (Particles[temp]->IsFree && (Particles[temp]->x(0) >= BC.OutFlowLoc) )
					{
						BC.OutPart.Push(temp);
						Particles[temp]->InOut = 2;
					}
					temp = Particles[temp]->LL;
				}
			}
		}
	}

	Vec3_t vel;
	double den;
	if (BC.InPart.Size()>0)
		#pragma omp parallel for schedule(static) private(vel,den) num_threads(Nproc)
		for (int i=0 ; i<BC.InPart.Size() ; i++)
		{
			size_t a = BC.InPart[i];
			InCon(Particles[a]->x,vel,den,BC);
			Particles[a]->v  = vel;
			Particles[a]->vb = vel;
			Particles[a]->va = vel;
			Particles[a]->Density  = den;
			Particles[a]->Densityb = den;
			Particles[a]->Densitya = den;
    			Particles[a]->Pressure = EOS(Particles[a]->PresEq, Particles[a]->Cs, Particles[a]->P0,Particles[a]->Density, Particles[a]->RefDensity);
		}

	double temp11;
	if (BC.MassConservation)
		temp11 = BC.InPart.Size()*1.0/(BC.OutPart.Size()*1.0);
	else
		temp11 = 1.0;

	if (BC.OutPart.Size()>0)
		#pragma omp parallel for schedule(static) private(vel,den) num_threads(Nproc)
		for (int i=0 ; i<BC.OutPart.Size() ; i++)
		{
			size_t a = BC.OutPart[i];
			OutCon(Particles[a]->x,vel,den,BC);
			if (norm(BC.outv)>0.0 || temp11 != 1.0)
			{
				Particles[a]->v  = temp11*vel;
				Particles[a]->vb = temp11*vel;
				Particles[a]->va = temp11*vel;
			}
			if (BC.outDensity>0.0)
			{
				Particles[a]->Density  = den;
				Particles[a]->Densityb = den;
 				Particles[a]->Densitya = den;
    				Particles[a]->Pressure = EOS(Particles[a]->PresEq, Particles[a]->Cs, Particles[a]->P0,Particles[a]->Density, Particles[a]->RefDensity);
			}
		}
}

inline void Domain::InFlowBCLeave()
{
	size_t a,b;
	Array <int> DelPart,TempPart;
	Array<std::pair<Vec3_t,size_t> > AddPart;

	#pragma omp parallel for schedule(static) num_threads(Nproc)
	for (int i=0; i<Particles.Size(); i++)
		if ((Particles[i]->x(0) > TRPR(0)) || (Particles[i]->x(1) > TRPR(1)) || (Particles[i]->x(2) > TRPR(2)) ||
				(Particles[i]->x(0) < BLPF(0)) || (Particles[i]->x(1) < BLPF(1)) || (Particles[i]->x(2) < BLPF(2)))
		{
			Particles[i]->InOut	= 0;
			omp_set_lock(&dom_lock);
			DelPart.Push(i);
			omp_unset_lock(&dom_lock);
		}

	if (BC.InOutFlow==1 || BC.InOutFlow==3)
	{
		for (size_t i=0 ; i<BC.InPart.Size() ; i++)
			if(Particles[BC.InPart[i]]->x(0) > BC.InFlowLoc1)
			{
				Vec3_t temp1	 = Particles[BC.InPart[i]]->x;
				temp1(0)	-= (BC.InFlowLoc3-BC.InFlowLoc2+InitialDist);
				Particles[BC.InPart[i]]->InOut		= 0;
				Particles[BC.InPart[i]]->FirstStep	= false;
				Particles[BC.InPart[i]]->ShepardCounter	= 0;
				AddPart.Push(std::make_pair(temp1,BC.InPart[i]));
				TempPart.Push(i);
			}
		BC.InPart.DelItems(TempPart);
		TempPart.Clear();
	}

	if (AddPart.Size() >= DelPart.Size())
	{
		for (size_t i=0 ; i<DelPart.Size() ; i++)
		{
			a = DelPart[i];
			b = AddPart[i].second;
			Particles[a]->x 		= AddPart[i].first;
			Particles[a]->Material		= 1;
			Particles[a]->InOut		= 1;
			Particles[a]->FirstStep		= false;

			Particles[a]->P0		= Particles[b]->P0;
			Particles[a]->PresEq		= Particles[b]->PresEq;
			Particles[a]->Cs		= Particles[b]->Cs;

			Particles[a]->Alpha		= Particles[b]->Alpha;
			Particles[a]->Beta		= Particles[b]->Beta;
			Particles[a]->Mu		= Particles[b]->Mu;
			Particles[a]->MuRef		= Particles[b]->MuRef;
			Particles[a]->T0		= 0.0; // Inflow is not capable of injecting non-Newtonian fluid

			Particles[a]->Mass 		= Particles[b]->Mass;
			Particles[a]->h			= Particles[b]->h;

			Particles[a]->ID 		= Particles[b]->ID;

			Particles[a]->TI		= 0.0; // Inflow is not capable of condidering the tensile instability

			Particles[a]->RefDensity	= Particles[b]->RefDensity; // The density for inflow must always be defined

			Particles[a]->ct		= Particles[b]->ct;

			Particles[a]->Shepard		= Particles[b]->Shepard;
			Particles[a]->ShepardStep	= Particles[b]->ShepardStep;
			Particles[a]->ShepardCounter	= Particles[b]->ShepardCounter;

			Particles[a]->LES		= Particles[b]->LES;
			Particles[a]->CSmag		= Particles[b]->CSmag;
			BC.InPart.Push(a);
		}

		if (AddPart.Size() != DelPart.Size())
			for (size_t i=DelPart.Size() ; i<AddPart.Size() ; i++)
			{
				b = AddPart[i].second;

				Particles.Push(new Particle(Particles[b]->ID,AddPart[i].first,Particles[b]->v,Particles[b]->Mass,Particles[b]->RefDensity,Particles[b]->h,false));

				a = Particles.Size()-1;
				Particles[a]->Material		= 1;
				Particles[a]->InOut		= 1;
				Particles[a]->FirstStep		= false;

				Particles[a]->P0		= Particles[b]->P0;
				Particles[a]->PresEq		= Particles[b]->PresEq;
				Particles[a]->Cs		= Particles[b]->Cs;

				Particles[a]->Alpha		= Particles[b]->Alpha;
				Particles[a]->Beta		= Particles[b]->Beta;
				Particles[a]->Mu		= Particles[b]->Mu;
				Particles[a]->MuRef		= Particles[b]->MuRef;
				Particles[a]->T0		= 0.0; // Inflow is not capable of injecting non-Newtonian fluid

				Particles[a]->TI		= 0.0; // Inflow is not capable of condidering the tensile instability

				Particles[a]->ct		= Particles[b]->ct;

				Particles[a]->Shepard		= Particles[b]->Shepard;
				Particles[a]->ShepardStep	= Particles[b]->ShepardStep;
				Particles[a]->ShepardCounter	= Particles[b]->ShepardCounter;

				Particles[a]->LES		= Particles[b]->LES;
				Particles[a]->CSmag		= Particles[b]->CSmag;
				BC.InPart.Push(a);
			}

		DelPart.Clear();
		AddPart.Clear();
	}
	else
	{
		for (size_t i=0 ; i<AddPart.Size() ; i++)
		{
			a = DelPart[i];
			b = AddPart[i].second;
			Particles[a]->x 		= AddPart[i].first;
			Particles[a]->Material		= 1;
			Particles[a]->InOut		= 1;
			Particles[a]->FirstStep		= false;

			Particles[a]->P0		= Particles[b]->P0;
			Particles[a]->PresEq		= Particles[b]->PresEq;
			Particles[a]->Cs		= Particles[b]->Cs;

			Particles[a]->Alpha		= Particles[b]->Alpha;
			Particles[a]->Beta		= Particles[b]->Beta;
			Particles[a]->Mu		= Particles[b]->Mu;
			Particles[a]->MuRef		= Particles[b]->MuRef;
			Particles[a]->T0		= 0.0; // Inflow is not capable of injecting non-Newtonian fluid

			Particles[a]->Mass 		= Particles[b]->Mass;
			Particles[a]->h			= Particles[b]->h;

			Particles[a]->ID 		= Particles[b]->ID;

			Particles[a]->TI		= 0.0; // Inflow is not capable of condidering the tensile instability

			Particles[a]->RefDensity	= Particles[b]->RefDensity; // The density for inflow must always be defined

			Particles[a]->ct		= Particles[b]->ct;

			Particles[a]->Shepard		= Particles[b]->Shepard;
			Particles[a]->ShepardStep	= Particles[b]->ShepardStep;
			Particles[a]->ShepardCounter	= Particles[b]->ShepardCounter;

			Particles[a]->LES		= Particles[b]->LES;
			Particles[a]->CSmag		= Particles[b]->CSmag;
			BC.InPart.Push(a);
		}
		for (size_t i=AddPart.Size() ; i<DelPart.Size() ; i++)
		{
			TempPart.Push(DelPart[i]);
		}
		Particles.DelItems(TempPart);
		BC.inoutcounter = 1;
		DelPart.Clear();
		AddPart.Clear();
	}
}

inline void Domain::CheckParticleLeave ()
{
	Array <int> DelParticles;

	#pragma omp parallel for schedule(static) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t i=0; i<Particles.Size(); i++){	//Like in Domain::Move
	#else
	for (int i=0; i<Particles.Size(); i++){//Like in Domain::Move
	#endif
    
		if ((Particles[i]->x(0) > TRPR(0)) || (Particles[i]->x(1) > TRPR(1)) || (Particles[i]->x(2) > TRPR(2)) ||
				(Particles[i]->x(0) < BLPF(0)) || (Particles[i]->x(1) < BLPF(1)) || (Particles[i]->x(2) < BLPF(2)))
		{
			omp_set_lock(&dom_lock);
			DelParticles.Push(i);
			omp_unset_lock(&dom_lock);
		}
    }

	if (DelParticles.Size()>0)
	{
		std::cout<< DelParticles.Size()<< " particle(s) left the Domain"<<std::endl;
		for (size_t i=0; i<DelParticles.Size(); i++)
		{
			std::cout<<""<<std::endl;
			std::cout<<"Particle Number   = "<<DelParticles[i]<<std::endl;
			std::cout<<"Particle Material = "<<Particles[DelParticles[i]]->Material<<std::endl;
			std::cout<<"x = "<<Particles[DelParticles[i]]->x<<std::endl;
			std::cout<<"v = "<<Particles[DelParticles[i]]->v<<std::endl;
			std::cout<<"a = "<<Particles[DelParticles[i]]->a<<std::endl;
		}
		Particles.DelItems(DelParticles);
	}
}

}; // namespace SPH