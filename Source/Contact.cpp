#include "matvec.h" 

namespace SPH {

inline Vec3_t cross(const Vec3_t &a, const Vec3_t &b){
	Vec3_t ret;
	ret(0) = a(1)*b(2)-a(2)*b(1);
	ret(1) = a(2)*b(0)-a(0)*b(2);
	ret(2) = a(0)*b(1)-a(1)*b(0);
	return ret;
}
	
void Domain::AddTrimeshParticles(const TriMesh &mesh, const float &hfac, const int &id){
	
	first_fem_particle_idx = Particles.Size();
	double Density =0.;
	double h;
	bool Fixed = true;	//Always are fixed ...
	contact_surf_id = id;
	trimesh = &mesh;
	
	for ( int e = 0; e < mesh.element.Size(); e++ ){
		Vec3_t pos = mesh.element[e]->centroid;
		h = hfac * mesh.element[e]->radius;
		Particles.Push(new Particle(id,pos,Vec3_t(0,0,0),0.0,Density,h,Fixed));
		Particles[first_fem_particle_idx + e] -> normal  = mesh.element[e] -> normal;
		Particles[first_fem_particle_idx + e] -> element = e; 
	}
	cout << Particles.Size() - first_fem_particle_idx << "particles added with ID " << contact_surf_id <<endl;
	cout << first_fem_particle_idx << " is the first solid particle index."<<endl;
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
			// cout << "Fist Particle index, ID: "<<P1 << ", " << Particles[P1]->ID  <<endl;
			// cout << "Sec Particle index, ID: "<<P2 << ", " << Particles[P2]->ID  <<endl;
			if (Particles[P1]->ID == contact_surf_id || Particles[P2]->ID == contact_surf_id ) {
				//cout << "contact surface id found"<<endl;
				if (Particles[P1]->ID == id_free_surf || Particles[P2]->ID == id_free_surf ) {
					Vec3_t xij	= Particles[P1]->x - Particles[P2]->x;
					double r = norm(xij);
					double rcutoff = ( Particles[P1]->h + Particles[P2]->h ) / 2.;
					//cout << "r, rcutoff, h1, h2"<< r << ", "<< rcutoff << ", "<< Particles[P1]->h <<", "<<Particles[P2]->h<<endl;
					if ( r < 2.0 *rcutoff ){
					//cout << "Found contact pair: "<< P1 << ", " << P2 << endl;
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
void Domain::CalcContactForces(){
	max_contact_force = 0.;
	int inside_pairs = 0;
	#pragma omp parallel for schedule (static) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t k=0; k<Nproc;k++) 
	#else
	for (int k=0; k<Nproc;k++) 
	#endif	
	{
		int P1,P2;
		Vec3_t xij;
		double h,K;
		// Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
		//IT IS CONVENIENT TO FIX SINCE FSMPairs are significantly smaller
		//cout << "Contact pair size: "<<ContPairs[k].Size()<<endl;
		for (size_t a = 0; a < ContPairs[k].Size();a++) {
			//P1 is SPH particle, P2 is CONTACT SURFACE (FEM) Particle
			if (Particles[ContPairs[k][a].first]->ID == contact_surf_id ) 	{ 	//Cont Surf is partcicles from FEM
				P1 = ContPairs[k][a].second; P2 = ContPairs[k][a].first; 	}
			else {
				P1 = ContPairs[k][a].first; P2 = ContPairs[k][a].second; } 
		
			Vec3_t vr = Particles[P1]->v - Particles[P2]->v;		//Fraser 3-137
			// cout << "Particle P1v: "<<Particles[P1]->v<<endl;
			// cout << "Particle P2v: "<<Particles[P2]->v<<endl;
			//ok, FEM Particles normals can be calculated by two ways, the one used to
			//calculate SPH ones, and could be given by mesh input
			//delta_ Is the projection of relative velocity 
			double delta_ = - dot( Particles[P2]->normal , vr);	//Penetration rate, Fraser 3-138
			
			// cout << "Particle Normal: "<<Particles[P2]->normal<<endl;
			// cout << " vr: "<< vr<<endl;
			// cout << "delta_: "<<delta_<<endl;
			//Check if SPH and fem particles are approaching each other
			if (delta_ > 0 ){
				Element* e = trimesh-> element[Particles[P2]->element];
				double pplane = e -> pplane; 
				double deltat_cont = ( Particles[P1]->h + pplane - dot (Particles[P2]->normal,	Particles[P1]->x) ) / (-delta_);								//Eq 3-142 
				Vec3_t Ri = Particles[P1]->x + deltat_cont * vr;	//Eq 3-139 Ray from SPH particle in the rel velocity direction
				//cout << "dt contact: "<<deltat_cont<<endl;
				//Check for contact in this time step 
				//Calculate time step for external forces
				
				if (deltat_cont < deltat){
					//cout << "Inside dt contact" <<endl;
					//Find point of contact Qj
					Vec3_t Qj = Particles[P1]->x + (Particles[P1]->v * deltat_cont) - ( Particles[P1]->h, Particles[P2]->normal); //Fraser 3-146
					//Check if it is inside triangular element
					//Find a vector 
					//Fraser 3-147
					bool inside = true;
					int i=0,j;			
					while (i<3 && inside){
						j = i+1;	if (j>2) j = 0;
						double crit = dot (cross ( *trimesh->node[e -> node[j]] - *trimesh->node[e -> node[i]],
																																Qj  - *trimesh->node[e -> node[i]]),
															Particles[P2]->normal);
						if (crit < 0.0) inside = false;
						i++;
					}
					
					if (inside ) { //Contact point inside element, contact proceeds
						//Recalculate vr (for large FEM mesh densities)
						//cout << "inside element"<<endl;
						
						//Calculate penetration depth
						double delta = (deltat - deltat_cont) * delta_;
						
						// DAMPING
						//Calculate SPH and FEM elements stiffness (series)
						//Since FEM is assumed as rigid, stiffness is simply the SPH one 
						double kij = PFAC * Particles[P1]-> cont_stiff;
						double omega = sqrt (kij/Particles[P1]->Mass);
						double psi_cont = 2. * Particles[P1]->Mass * omega * DFAC; // Fraser Eqn 3-158
						
						// TANGENTIAL COMPONENNT
						//Fraser Eqn 3-167
						
						double force2 = dot(Particles[P1] -> contforce,Particles[P1] -> contforce);
						omp_set_lock(&Particles[P1]->my_lock);
						Particles[P1] -> contforce = (kij * delta - psi_cont * delta_) * Particles[P2]->normal; // NORMAL DIRECTION
						Particles[P1] -> a += Particles[P1] -> contforce / Particles[P1] -> Mass; 
						omp_unset_lock(&Particles[P1]->my_lock);
						
						if (force2 > max_contact_force) max_contact_force = force2;
						inside_pairs++;
					}// if inside
				} //deltat <min

			}//delta_ > 0 : PARTICLES ARE APPROACHING EACH OTHER
		}//Contact Pairs
	}//Nproc
	max_contact_force = sqrt (max_contact_force);
	cout << "Max Contact Force: "<< max_contact_force << "Time: " << Time << ", Pairs"<<inside_pairs<<endl;
}

}; //SPH