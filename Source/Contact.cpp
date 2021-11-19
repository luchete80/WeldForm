#include "matvec.h" 

namespace SPH {
	
void Domain::AddTrimeshParticles(const TriMesh &mesh, const int &id){
	
	double Density =0.;
	double h = 1.1;
	bool Fixed = true;	//Always are fixed ...
	contact_surf_id = id;
	 
	for ( int e = 0; e < mesh.element.Size(); e++ ){
		Vec3_t pos = mesh.element[e]->centroid;
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
			if (Particles[P1]->ID == contact_surf_id || Particles[P1]->ID == contact_surf_id ) {
				ContPairs[k].Push(std::make_pair(P1, P2));
				//ContPairs[k].Push(FSMPairs[k][a]);
				//If the problem is not thermal (only mechanic)
				//Could be deleted this pair in Whole Pairs
				
				FSMPairs[k].DelItem(a);//Erase, NOT EFFICIENT
				a--;
			}
		}
	
	}
}

}; //SPH