#include "Domain.h"

#include <CompactNSearch> //NEW WAY

using namespace std;
namespace SPH {
	
inline void Domain::CellInitiate () {
	if (!(norm(TRPR)>0.0) && !(norm(BLPF)>0.0))
	{
		// Calculate Domain Size
		BLPF = Particles[0]->x;
		TRPR = Particles[0]->x;
		hmax = Particles[0]->h;
		rhomax = Particles[0]->Density;

		for (size_t i=0; i<Particles.Size(); i++)
		{
			if (Particles[i]->x(0) > TRPR(0)) TRPR(0) = Particles[i]->x(0);
			if (Particles[i]->x(1) > TRPR(1)) TRPR(1) = Particles[i]->x(1);
			if (Particles[i]->x(2) > TRPR(2)) TRPR(2) = Particles[i]->x(2);

			if (Particles[i]->x(0) < BLPF(0)) BLPF(0) = Particles[i]->x(0);
			if (Particles[i]->x(1) < BLPF(1)) BLPF(1) = Particles[i]->x(1);
			if (Particles[i]->x(2) < BLPF(2)) BLPF(2) = Particles[i]->x(2);

			if (Particles[i]->h > hmax) hmax=Particles[i]->h;
			if (Particles[i]->Density > rhomax) rhomax=Particles[i]->Density;
			if (Particles[i]->Mu > MuMax) MuMax=Particles[i]->Mu;
			if (Particles[i]->Cs > CsMax) CsMax=Particles[i]->Cs;
		}
	}

	// Override the calculated domain size
	if (DomMax(0)>TRPR(0)) TRPR(0) = DomMax(0);
	if (DomMax(1)>TRPR(1)) TRPR(1) = DomMax(1);
	if (DomMax(2)>TRPR(2)) TRPR(2) = DomMax(2);
	if (DomMin(0)<BLPF(0)) BLPF(0) = DomMin(0);
	if (DomMin(1)<BLPF(1)) BLPF(1) = DomMin(1);
	if (DomMin(2)<BLPF(2)) BLPF(2) = DomMin(2);


	//Because of Hexagonal close packing in x direction domain is modified
	if (!BC.Periodic[0]) {TRPR(0) += hmax/2;	BLPF(0) -= hmax/2;}else{TRPR(0) += R; BLPF(0) -= R;}
	if (!BC.Periodic[1]) {TRPR(1) += hmax/2;	BLPF(1) -= hmax/2;}else{TRPR(1) += R; BLPF(1) -= R;}
	if (!BC.Periodic[2]) {TRPR(2) += hmax/2;	BLPF(2) -= hmax/2;}else{TRPR(2) += R; BLPF(2) -= R;}

    // Calculate Cells Properties
	switch (Dimension)
	{case 2:
		if (double (ceil(((TRPR(0)-BLPF(0))/(Cellfac*hmax)))-((TRPR(0)-BLPF(0))/(Cellfac*hmax)))<(hmax/10.0))
			CellNo[0] = int(ceil((TRPR(0)-BLPF(0))/(Cellfac*hmax)));
		else
			CellNo[0] = int(floor((TRPR(0)-BLPF(0))/(Cellfac*hmax)));

		if (double (ceil(((TRPR(1)-BLPF(1))/(Cellfac*hmax)))-((TRPR(1)-BLPF(1))/(Cellfac*hmax)))<(hmax/10.0))
			CellNo[1] = int(ceil((TRPR(1)-BLPF(1))/(Cellfac*hmax)));
		else
			CellNo[1] = int(floor((TRPR(1)-BLPF(1))/(Cellfac*hmax)));

		CellNo[2] = 1;

		CellSize  = Vec3_t ((TRPR(0)-BLPF(0))/CellNo[0],(TRPR(1)-BLPF(1))/CellNo[1],0.0);
		break;

	case 3:
		if (double (ceil(((TRPR(0)-BLPF(0))/(Cellfac*hmax)))-((TRPR(0)-BLPF(0))/(Cellfac*hmax)))<(hmax/10.0))
			CellNo[0] = int(ceil((TRPR(0)-BLPF(0))/(Cellfac*hmax)));
		else
			CellNo[0] = int(floor((TRPR(0)-BLPF(0))/(Cellfac*hmax)));

		if (double (ceil(((TRPR(1)-BLPF(1))/(Cellfac*hmax)))-((TRPR(1)-BLPF(1))/(Cellfac*hmax)))<(hmax/10.0))
			CellNo[1] = int(ceil((TRPR(1)-BLPF(1))/(Cellfac*hmax)));
		else
			CellNo[1] = int(floor((TRPR(1)-BLPF(1))/(Cellfac*hmax)));

		if (double (ceil(((TRPR(2)-BLPF(2))/(Cellfac*hmax)))-((TRPR(2)-BLPF(2))/(Cellfac*hmax)))<(hmax/10.0))
			CellNo[2] = int(ceil((TRPR(2)-BLPF(2))/(Cellfac*hmax)));
		else
			CellNo[2] = int(floor((TRPR(2)-BLPF(2))/(Cellfac*hmax)));

		CellSize  = Vec3_t ((TRPR(0)-BLPF(0))/CellNo[0],(TRPR(1)-BLPF(1))/CellNo[1],(TRPR(2)-BLPF(2))/CellNo[2]);
		break;

	default:
    	std::cout << "Please correct the dimension (2=>2D or 3=>3D) and run again" << std::endl;
		abort();
		break;
	}

	// Periodic BC modifications
	if (BC.Periodic[0]) CellNo[0] += 2;
    if (BC.Periodic[1]) CellNo[1] += 2;
    if (BC.Periodic[2]) CellNo[2] += 2;

    if (BC.Periodic[0]) DomSize[0] = (TRPR(0)-BLPF(0));
    if (BC.Periodic[1]) DomSize[1] = (TRPR(1)-BLPF(1));
    if (BC.Periodic[2]) DomSize[2] = (TRPR(2)-BLPF(2));

    // Initiate Head of Chain array for Linked-List
    HOC = new int**[(int) CellNo[0]];
    for(int i =0; i<CellNo[0]; i++){
       HOC[i] = new int*[CellNo[1]];
       for(int j =0; j<CellNo[1]; j++){
           HOC[i][j] = new int[CellNo[2]];
           for(int k = 0; k<CellNo[2];k++){
              HOC[i][j][k] = -1;
           }
       }
    }
    // Initiate Pairs array for neibour searching
    for(size_t i=0 ; i<Nproc ; i++) {
			SMPairs.Push(Initial);
			NSMPairs.Push(Initial);
			FSMPairs.Push(Initial);
			RIGPairs.Push(Initial);
			
			ContPairs.Push(Initial);
      
      //New integration/sum
      Array <size_t> a;
      ilist_SM.Push(a);jlist_SM.Push(a);
      //ipair_SM.Push(a);jpair_SM.Push(a);
    }
}

inline void Domain::ListGenerate ()
{
	int i, j, k, temp=0;
	switch (Dimension)
	{case 2:
		for (size_t a=0; a<Particles.Size(); a++)
		{
			i= (int) (floor((Particles[a]->x(0) - BLPF(0)) / CellSize(0)));
			j= (int) (floor((Particles[a]->x(1) - BLPF(1)) / CellSize(1)));

			if (i<0)
            {
                    if ((BLPF(0) - Particles[a]->x(0)) <= hmax) i=0;
                            else std::cout<<"Leaving i<0"<<std::endl;
            }
            if (j<0)
            {
                    if ((BLPF(1) - Particles[a]->x(1)) <= hmax) j=0;
                            else std::cout<<"Leaving j<0"<<std::endl;
            }
			if (i>=CellNo[0])
			{
					if ((Particles[a]->x(0) - TRPR(0)) <= hmax) i=CellNo[0]-1;
							else std::cout<<"Leaving i>=CellNo"<<std::endl;
			}
            if (j>=CellNo[1])
            {
                    if ((Particles[a]->x(1) - TRPR(1)) <= hmax) j=CellNo[1]-1;
                            else std::cout<<"Leaving j>=CellNo"<<std::endl;
            }

			temp = HOC[i][j][0];
			HOC[i][j][0] = a;
			Particles[a]->LL = temp;
			Particles[a]->CC[0] = i;
			Particles[a]->CC[1] = j;
			Particles[a]->CC[2] = 0;
			if (!Particles[a]->IsFree) FixedParticles.Push(a);
		}
		break;

	case 3:
		for (size_t a=0; a<Particles.Size(); a++)
		{
			i= (int) (floor((Particles[a]->x(0) - BLPF(0)) / CellSize(0)));
			j= (int) (floor((Particles[a]->x(1) - BLPF(1)) / CellSize(1)));
			k= (int) (floor((Particles[a]->x(2) - BLPF(2)) / CellSize(2)));

            if (i<0) i = 0;
            // {
                    // if ((BLPF(0) - Particles[a]->x(0))<=hmax) i=0;
                            // else std::cout<<"Leaving, particle "<<a<< "yield "<<Particles[a]->Sigmay<<", eff_str_rate "<<Particles[a]->eff_strain_rate<<std::endl;
            // }
            if (j<0) j = 0 ;
            // {
                    // if ((BLPF(1) - Particles[a]->x(1))<=hmax) j=0;
                            // else std::cout<<"Leaving particle "<<a<<"yield "<<Particles[a]->Sigmay<<", eff_str_rate "<<Particles[a]->eff_strain_rate<<std::endl;
            // }
            if (k<0) k = 0;
            // {
                    // if ((BLPF(2) - Particles[a]->x(2))<=hmax) k=0;
                            // else std::cout<<"Leaving particle"<< a << "yield "<<Particles[a]->Sigmay<<", eff_str_rate "<<Particles[a]->eff_strain_rate<<std::endl;
            // }
			if (i>=CellNo[0])i=CellNo[0]-1;
			// {
					// if ((Particles[a]->x(0) - TRPR(0))<=hmax) i=CellNo[0]-1;
							// else std::cout<<"Leavin particle "<<a<< "yield "<<Particles[a]->Sigmay<<", eff_str_rate "<<Particles[a]->eff_strain_rate<<std::endl;
			// }
            if (j>=CellNo[1])j=CellNo[1]-1;
            // {
                    // if ((Particles[a]->x(1) - TRPR(1))<=hmax) j=CellNo[1]-1;
                            // else std::cout<<"Leaving particle "<<a<<"yield "<<Particles[a]->Sigmay<<", eff_str_rate "<<Particles[a]->eff_strain_rate<<std::endl;
            // }
            if (k>=CellNo[2])k=CellNo[2]-1;
            // {
                    // if ((Particles[a]->x(2) - TRPR(2))<=hmax) k=CellNo[2]-1;
                            // else std::cout<<"Leaving particle"<<a<<"yield "<<Particles[a]->Sigmay<<", eff_str_rate "<<Particles[a]->eff_strain_rate<<std::endl;
            // }

            temp = HOC[i][j][k];
			HOC[i][j][k] = a;
			Particles[a]->LL = temp;
			Particles[a]->CC[0] = i;
			Particles[a]->CC[1] = j;
			Particles[a]->CC[2] = k;
			if (!Particles[a]->IsFree) FixedParticles.Push(a);
		}
		break;

	default:
    	std::cout << "Please correct the dimension (2=>2D or 3=>3D) and run again" << std::endl;
		abort();
		break;
	}

	if (BC.Periodic[0]) {
	   for(int j =0; j<CellNo[1]; j++)
		   for(int k =0; k<CellNo[2]; k++) {
			  HOC[CellNo[0]-1][j][k] =  HOC[1][j][k];
			  HOC[CellNo[0]-2][j][k] =  HOC[0][j][k];
		   }
	} 
	if (BC.Periodic[1]) {
	   for(int i =0; i<CellNo[0]; i++)
		   for(int k =0; k<CellNo[2]; k++) {
			  HOC[i][CellNo[1]-1][k] =  HOC[i][1][k];
			  HOC[i][CellNo[1]-2][k] =  HOC[i][0][k];
		   }
	}
	if (BC.Periodic[2]) {
	   for(int i =0; i<CellNo[0]; i++)
		   for(int j =0; j<CellNo[1]; j++) {
				  HOC[i][j][CellNo[2]-1] =  HOC[i][j][1];
				  HOC[i][j][CellNo[2]-2] =  HOC[i][j][0];
			   }
	}
}

inline void Domain::CellReset ()
{

    #pragma omp parallel for schedule (static) num_threads(Nproc)

    for(int i =0; i<CellNo[0]; i++)
    {
		for(int j =0; j<CellNo[1]; j++)
		for(int k =0; k<CellNo[2];k++)
		{
			HOC[i][j][k] = -1;
		}
    }
	#pragma omp parallel for schedule(static) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t a=0; a<Particles.Size(); a++)	//Like in Domain::Move
	#else
	for (int a=0; a<Particles.Size(); a++)//Like in Domain::Move
	#endif
	{

    	Particles[a]->LL = -1;
    }

    FixedParticles.Clear();
}

using namespace CompactNSearch;

inline void Domain::MainNeighbourSearch_CNS(const double &radius){
  
  std::vector<std::array<Real, 3>> positions;


	for (unsigned int p = 0; p < Particles.Size(); p++){
		
		std::array<Real, 3> x = {{
			static_cast<Real>(Particles[p]->x[0]),
			static_cast<Real>(Particles[p]->x[1]),
			static_cast<Real>(Particles[p]->x[2])
		}};
		positions.push_back(x);

	}
		
	//std::random_shuffle(positions.begin(), positions.end());
	cout << "Finding NEighbours"<<endl;
	NeighborhoodSearch nsearch(radius, true);
	nsearch.add_point_set(positions.front().data(), positions.size(), true, true);
	nsearch.add_point_set(positions.front().data(), positions.size(), true, true);
	nsearch.find_neighbors();

	nsearch.update_point_sets();
	std::vector<std::vector<unsigned int>> neighbors2;
	nsearch.find_neighbors(0, 1, neighbors2);
	std::vector<std::vector<unsigned int>> neighbors3;
	nsearch.find_neighbors(1, 2, neighbors3);


	ofstream outfind2; // outdata is like cin
	outfind2.open("find2_part.txt"); // opens the file
	//Pass to domain
	std::set< std:: pair<int,int> > neigbours_set;
	auto const& d = nsearch.point_set(0);
	for (int i = 0; i < d.n_points(); ++i){
		const std::vector<unsigned int>& nbs = d.neighbor_list(0, i);
		//res += static_cast<unsigned long>(d.n_neighbors(0, i));		

		for (int k=0;k< nbs.size();k++) {
			outfind2<< i << ", "<<nbs[k]<<endl;
			neigbours_set.insert(std::make_pair( std::min(i,int(nbs[k])), std::max(i,int(nbs[k]))) );
			
		}
	}
  
    //TEST
	// nsearch.z_sort();
	// for (auto i = 0u; i < nsearch.n_point_sets(); ++i)
	// {
		// auto const& d = nsearch.point_set(i);
		// d.sort_field(positions.data());

	// }
	// nsearch.find_neighbors();
	
}

inline void Domain::MainNeighbourSearch() {
    int q1;
		//cout << "id free surf"<<id_free_surf<<endl;
    if (BC.Periodic[0]) {
	#pragma omp parallel for schedule (dynamic) num_threads(Nproc)
	for (q1=1;q1<(CellNo[0]-1); q1++)	YZPlaneCellsNeighbourSearch(q1);
    } else {
	#pragma omp parallel for schedule (dynamic) num_threads(Nproc)
    	for (q1=0;q1<CellNo[0]; q1++)	YZPlaneCellsNeighbourSearch(q1);
    }
	m_isNbDataCleared = false;
}

inline bool  Domain::CheckRadius(Particle* P1, Particle *P2){
	bool ret = false;
	
	double h	= (P1->h+P2->h)/2;
	Vec3_t xij	= P1->x - P2->x;
	Periodic_X_Correction(xij, h, P1, P2);
	double rij	= norm(xij);
	if ((rij/h)<=Cellfac) ret = true;
	// cout << "Checking radius "<<endl;
	// cout << "rij h rij/h cellfac "<<rij<<", "<< h << ", " << rij/h<<", "<<Cellfac<<endl;
	return ret;
}

inline void Domain::AllocateNbPair(const int &temp1, const int &temp2, const int &T){

		if (contact) { //ONLY PAIRS WITH ONE PARTICLE ON THE SURFACE AND ONE PARTICLE ON THE RIGID SURFACE
			// if (  (Particles[temp1]->ID == contact_surf_id && Particles[temp2]->ID == id_free_surf ) || 
              // (Particles[temp1]->ID == id_free_surf && Particles[temp2]->ID == contact_surf_id )) {
      bool is_contact = false;
      for (int m=0;m<meshcount;m++){
        if (Particles[temp1]->ID == contact_surf_id[m] || Particles[temp2]->ID == contact_surf_id[m] )
          is_contact = true;
      }
      if (is_contact)  {
        if (Particles[temp1]->ID == id_free_surf || Particles[temp2]->ID == id_free_surf ) {
            //ORIGINAL OLD
            RIGPairs[T].Push(std::make_pair(temp1, temp2));
            
            //NEW
            // if ( dot(Particles[temp1]->x - Particles[temp2]->x,Particles[temp1]->x - Particles[temp2]->x) < 
                  // ( Particles[temp1]->h + Particles[temp2]->h ) *( Particles[temp1]->h + Particles[temp2]->h ))
              // ContPairs[T].Push(std::make_pair(temp1, temp2));

        }
        return;	//If either one particle or another is in the surface 
      } 
      
		}
		int i,j;
		if ( CheckRadius(Particles[temp1],Particles[temp2])){
			if (Particles[temp1]->IsFree || Particles[temp2]->IsFree) {
				if (Particles[temp1]->Material == Particles[temp2]->Material)
				{
					if (Particles[temp1]->IsFree*Particles[temp2]->IsFree) {//Both free, most common
						SMPairs[T].Push(std::make_pair(temp1, temp2));
            // #ifdef NONLOCK_SUM
            // i = std::min(temp1,temp2);
            // j = std::max(temp1,temp2);
            // ilist_SM[T].Push(i); //THIS COULD BE DONE AFTER 
            // jlist_SM[T].Push(j);
            // ipair_SM[T][i]++;
            // jpair_SM[T][j]++;
            // #endif
          }
					else
						FSMPairs[T].Push(std::make_pair(temp1, temp2)); //TEMPORARY
				} else
					NSMPairs[T].Push(std::make_pair(temp1, temp2));
			}
	}
}
inline void Domain::YZPlaneCellsNeighbourSearch(int q1) {
	int q3,q2;
	size_t T = omp_get_thread_num();

	for (BC.Periodic[2] ? q3=1 : q3=0;BC.Periodic[2] ? (q3<(CellNo[2]-1)) : (q3<CellNo[2]); q3++)
	for (BC.Periodic[1] ? q2=1 : q2=0;BC.Periodic[1] ? (q2<(CellNo[1]-1)) : (q2<CellNo[1]); q2++) {
		if (HOC[q1][q2][q3]==-1) continue;
		else {
			int temp1, temp2;
			temp1 = HOC[q1][q2][q3];

			while (temp1 != -1) {// The current cell  => self cell interactions
				temp2 = Particles[temp1]->LL;
				while (temp2 != -1){
						AllocateNbPair(temp1,temp2,T);
						temp2 = Particles[temp2]->LL;
				}//while

				// (q1 + 1, q2 , q3)
				if (q1+1< CellNo[0]) {
					temp2 = HOC[q1+1][q2][q3];
					while (temp2 != -1) {
						AllocateNbPair(temp1,temp2,T);
						temp2 = Particles[temp2]->LL;
					}//while temp2!=-1
				}// (q1 + 1, q2 , q3)

				// (q1 + a, q2 + 1, q3) & a[-1,1]
				if (q2+1< CellNo[1]) {
					for (int i = q1-1; i <= q1+1; i++) {
						if (i<CellNo[0] && i>=0) {
							temp2 = HOC[i][q2+1][q3];
							while (temp2 != -1)
							{
								AllocateNbPair(temp1,temp2,T);
								temp2 = Particles[temp2]->LL;
							}
						}
					}
				}

				// (q1 + a, q2 + b, q3 + 1) & a,b[-1,1] => all 9 cells above the current cell
				if (q3+1< CellNo[2]) {
					for (int j=q2-1; j<=q2+1; j++)
					for (int i=q1-1; i<=q1+1; i++) {
						if (i<CellNo[0] && i>=0 && j<CellNo[1] && j>=0) {
							temp2 = HOC[i][j][q3+1];
							while (temp2 != -1)
							{
								AllocateNbPair(temp1,temp2,T);
								temp2 = Particles[temp2]->LL;
							}
						}
					}
				}
				temp1 = Particles[temp1]->LL;
			}//while temp1 !=-1
		}
	}
}

inline void Domain::ClearNbData(){

	#pragma omp parallel for schedule (dynamic) num_threads(Nproc)
	for (int i=0 ; i<Nproc ; i++) { //In the original version this was calculated after
		SMPairs[i].Clear();
		FSMPairs[i].Clear();
		NSMPairs[i].Clear();
		RIGPairs[i].Clear();
		ContPairs[i].Clear();//New
	}
	CellReset();
	ListGenerate();
	m_isNbDataCleared = true;
}

inline void Domain::SaveNeighbourData(){
		std::vector <int> nb(Particles.Size());
		std::vector <int> contnb(Particles.Size());

	#pragma omp parallel for schedule (static) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t k=0; k<Nproc;k++) 
	#else
	for (int k=0; k<Nproc;k++) 
	#endif			
	{
			for (size_t a=0; a<SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			//cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
				nb[SMPairs[k][a].first ]+=1;
				nb[SMPairs[k][a].second]+=1;
				
			}
	}
		// for ( size_t k = 0; k < Nproc ; k++) {
			// for (size_t a=0; a<ContPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			// //cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
				// contnb[ContPairs[k][a].first ]+=1;
				// contnb[ContPairs[k][a].second]+=1;
				
			// }			
		// }

	#pragma omp parallel for schedule (static) num_threads(Nproc)		
	for (int p=0;p<Particles.Size();p++){
		omp_set_lock(&Particles[p]->my_lock);
		Particles[p]->Nb = nb[p];
		omp_unset_lock(&Particles[p]->my_lock);
	}
}

inline void Domain::SaveContNeighbourData(){
		std::vector <int> nb(Particles.Size());
		std::vector <int> contnb(Particles.Size());

		for ( size_t k = 0; k < Nproc ; k++) {
			for (size_t a=0; a<ContPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			//cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
				contnb[ContPairs[k][a].first ]+=1;
				contnb[ContPairs[k][a].second]+=1;
				
			}			
		}
		
		for (int p=0;p<Particles.Size();p++){
			if (p < first_fem_particle_idx[0])
				Particles[p]->ContNb = contnb[p];
		}
}

int Domain::AvgNeighbourCount(){	
		std::vector<int> nbcount(Particles.Size());
		
		#pragma omp parallel for schedule (static) num_threads(Nproc)
		for (int p=0;p<Nproc;p++)
			for (int i=0;i<SMPairs[p].size();i++){
				nbcount[SMPairs[p][i].first]++;
				nbcount[SMPairs[p][i].second]++;
			}
		int tot=0;
		for (int p=0;p<nbcount.size();p++)
		tot+=nbcount[p];
	
	return tot/Particles.Size();

}

// Calculate position ipl
inline void Domain::CalcPairPosList(){                             //Calculate position list for every particle
  
  int icount[Nproc],jcount[Nproc];
  for (int p=0;p<Nproc;p++) {
    icount[p]=jcount[p]=0;
  }
  
  #pragma omp parallel for schedule (static) num_threads(Nproc)
  for (int p=0;p<Nproc;p++){
    for (int i=0;i<Particles.Size();i++){
      ipl_SM[p][i] = icount[p];
      icount[p]    += ipair_SM[p][i];
      jpl_SM[p][i] = jcount[p];
      jcount[p]    += jpair_SM[p][i];
    }
  }
}

// Nishimura (2011)
// Algorithm 1: Creating contact candidate pair list and reference table from the table of neighbor particles
// 1: for i = 0 to N
// p  1 do {parallelized loop over particle index; Np is the total number of particles}
// 2: for Itable = 0 to njgi½i  1 do {loop over neighbor particle array}
// 3: Ilist  sjgi½i  1 + Itable {index of contact candidate pair list defined by Eq. (4)}
// 4: j  Anei½i; Itable
// 5: L
// p½Ilist  i {creating contact candidate pair list by Step 4}
// 6: Ln½Ilist  j
// 7: Aref½i; Itable  Ilist {Step (i) of Fig. 7}
// 8: for k = 0 to njli½j  1 do {Step (ii) of Fig. 7}
// 9: if Anei½j; N  1  k = i then {Step (iii) of Fig. 7}
// 10: Aref½j; njgi½j þ k  Ilist {Step (iv) of Fig. 7}
// 11: end if
// 12: end do
// 13: end do
// 14: end do

inline void Domain::CalcRefTable(){
  #pragma omp parallel for schedule (static) num_threads(Nproc)
  for (int p=0;p<Nproc;p++){
    for (int i = 0;i<Particles.Size();i++){
    size_t T = omp_get_thread_num();
      for (int n=0;n<ipl_SM[T][i];n++){
        //Aref = ipl_SM
      }//nb
    }//particle
  } //Thread
  
}

}; //SPH
