#include "matvec.h" 

namespace SPH {
	
void Domain::AddTrimeshParticles(const TriMesh &mesh, const float &hfac, const int &id){
	
	first_fem_particle_idx = Particles.Size();
	double Density =0.;
	double h;
	bool Fixed = true;	//Always are fixed ...
	contact_surf_id = id;
	 
	for ( int e = 0; e < mesh.element.Size(); e++ ){
		Vec3_t pos = mesh.element[e]->centroid;
		h = hfac * mesh.element[e]->radius;
		Particles.Push(new Particle(id,pos,Vec3_t(0,0,0),0.0,Density,h,Fixed));
	}
}

inline void Domain::ContactNbSearch(){
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
		//IT IS CONVENIENT TO FIX SINCE FSMPairs are significantly smaller
		for (size_t a=0; a<FSMPairs[k].Size();a++) {
			P1	= FSMPairs[k][a].first;
			P2	= FSMPairs[k][a].second;	
			if (Particles[P1]->ID == contact_surf_id || Particles[P2]->ID == contact_surf_id ) {
				if (Particles[P1]->ID == id_free_surf || Particles[P2]->ID == id_free_surf ) {
					Vec3_t xij	= Particles[P1]->x - Particles[P2]->x;
					double r = norm(xij);
					double rcutoff = ( Particles[P1]->h + Particles[P2]->h ) / 2.;
					//cout << "r, rcutoff, h1, h2"<< r << ", "<< rcutoff << ", "<< Particles[P1]->h <<", "<<Particles[P2]->h<<endl;
					if ( r < 2.0 *rcutoff ){
					cout << "Found contact pair: "<< P1 << ", " << P2 << endl;
					//ContPairs[k].Push(std::make_pair(P1, P2));
					ContPairs[k].Push(FSMPairs[k][a]);
					//If the problem is not thermal (only mechanic)
					//Could be deleted this pair in Whole Pairs
					}	//belongs to free surf
				}
				// Pair is removed either way is inside cutoff radius or is or in!
				FSMPairs[k].DelItem(a);//Erase, NOT EFFICIENT
				a--;
			}
		}

		// for (size_t a=0; a<SMPairs[k].Size();a++) {
			// P1	= SMPairs[k][a].first;
			// P2	= SMPairs[k][a].second;	
			// if (Particles[P1]->ID == contact_surf_id || Particles[P2]->ID == contact_surf_id ) {
				// cout << "Found contact pair: "<< P1 << ", " << P2 << endl;
				// //ContPairs[k].Push(std::make_pair(P1, P2));
				// ContPairs[k].Push(SMPairs[k][a]);
				// //If the problem is not thermal (only mechanic)
				// //Could be deleted this pair in Whole Pairs
				
				// SMPairs[k].DelItem(a);//Erase, NOT EFFICIENT
				// a--;
			// }
		// }
		
	}
}

//////////////////////////////// 
//// 
////////////////////////////////
void Domain::DetectContactPoints(){
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t k=0; k<Nproc;k++) 
	#else
	for (int k=0; k<Nproc;k++) 
	#endif	
	{
		Particle* P1,P2;
		Vec3_t xij;
		double h,K;
		// Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
		//IT IS CONVENIENT TO FIX SINCE FSMPairs are significantly smaller
		for (size_t a = 0; a < ContPairs[k].Size();a++) {
			P1	= Particles [ContPairs[k][a].first];
			P2	= Particles [ContPairs[k][a].second];				
			Particle* pi, pj; //SPH and contact part
			if (P1->ID == contact_surf_id ) { pi = P2; pj = P1; }
			else														{ pi = P1; pj = P2; }
		
			Vec3_t vr = pi->v - pj->v;		//Fraser 3-137
			double delta = - dot( , vr);	//Penetration rate, Fraser 3-138
		}
	}//Nproc
}

}; //SPH