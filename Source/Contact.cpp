#include "matvec.h" 
#define VMIN_FOR_FRICTION	1.e-3

namespace SPH {

	
void Domain::AddTrimeshParticles(const TriMesh &mesh, const float &hfac, const int &id){
	
	first_fem_particle_idx = Particles.Size();
	double Density =0.;
	double h;
	bool Fixed = false;	//Always are fixed ...
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
	//cout << "Performing Nb Search"<<endl;
	cont_pairs = 0;
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
		//cout << "Rig Pair size"<<RIGPairs[k].Size()<<endl;
		for (size_t a=0; a<RIGPairs[k].Size();a++) {
			
			P1	= RIGPairs[k][a].first;
			P2	= RIGPairs[k][a].second;	
			// cout << "Fist Particle index, ID: "<<P1 << ", " << Particles[P1]->ID  <<endl;
			// cout << "Sec Particle index, ID: "<<P2 << ", " << Particles[P2]->ID  <<endl;
			if (Particles[P1]->ID == contact_surf_id || Particles[P2]->ID == contact_surf_id ) {
				if (Particles[P1]->ID == id_free_surf || Particles[P2]->ID == id_free_surf ) {
					Vec3_t xij	= Particles[P1]->x - Particles[P2]->x;
					double r = norm(xij);
					double rcutoff = ( Particles[P1]->h + Particles[P2]->h ) / 2.;
					//cout << "r, rcutoff, h1, h2"<< r << ", "<< rcutoff << ", "<< Particles[P1]->h <<", "<<Particles[P2]->h<<endl;
					if ( r < 2.0 *rcutoff ){
					//cout << "Found contact pair: "<< P1 << ", " << P2 << endl;
					//ContPairs[k].Push(std::make_pair(P1, P2));
					ContPairs[k].Push(RIGPairs[k][a]);
					cont_pairs++;
					//If the problem is not thermal (only mechanic)
					//Could be deleted this pair in Whole Pairs
					}	//belongs to free surf
				}
				// Pair is removed either way is inside cutoff radius or is or in!
				RIGPairs[k].DelItem(a);//Erase, NOT EFFICIENT
				a--;
			}
		}
		//cout << "Found "<<cont_pairs<< " contact pairs."<<endl;
		
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
	// cout << "Contact Pairs Count"<<endl;
	// for (int k=0; k<Nproc;k++) 
		// cout << ContPairs[k].Size()<<", ";
	 // cout <<endl;
}

//////////////////////////////// 
//// 
////////////////////////////////
void Domain::CalcContactForces(){
	
	// #pragma omp parallel for num_threads(Nproc)
	// #ifdef __GNUC__
	// for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
	// #else
	// for (int i=0; i<Particles.Size(); i++)//Like in Domain::Move
	// #endif
	// {
		// Particles[i] -> contforce = 0.;
	// }
	double min_force_ts_=1000.;
// https://stackoverflow.com/questions/10850155/whats-the-difference-between-static-and-dynamic-schedule-in-openmp

	////////////////////////
	// DEBUG THINGS //////
	int inside_part[Particles.Size()];
	int inside_part_count = 0;
	double min_delta,max_delta;
	min_delta = 1000.; max_delta = 0.;
	int inside_time;
	
  //#pragma omp parallel for num_threads(Nproc)
  for (int i = 0;i<Particles.Size();i++){
		omp_set_lock(&Particles[i]->my_lock);
		Particles[i] -> contforce = 0.; //RESET    
		Particles[i] -> delta_cont = 0.; //RESET    
		Particles[i] -> tgdir = 0.;				//TODO: DELETE (DEBUG) 
		omp_unset_lock(&Particles[i]->my_lock);
		inside_part[i] = 0;
		inside_time=0;
  }
 
	
	max_contact_force = 0.;
	double min_contact_force = 1000.;
	int inside_pairs = 0;
	double force2 = 0.;
	double delta_ = 0.;
	double deltat_cont;
	
	int P1,P2;
	
	//#pragma omp parallel for schedule (static) private(force2,delta_,delta_cont, inside,i,j,crit) num_threads(Nproc)
	#ifdef __GNUC__
	for (size_t k=0; k<Nproc;k++) 
	#else
	for (int k=0; k<Nproc;k++) 
	#endif	
	{
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
			//cout << "Particle P1v: "<<Particles[P1]->v<<endl;
			//cout << "Particle P2v: "<<Particles[P2]->v<<endl;
			//ok, FEM Particles normals can be calculated by two ways, the one used to
			//calculate SPH ones, and could be given by mesh input
			//delta_ Is the projection of relative velocity 
			delta_ = - dot( Particles[P2]->normal , vr);	//Penetration rate, Fraser 3-138
			

			//Check if SPH and fem particles are approaching each other
			if (delta_ > 0 ){
				Element* e = trimesh-> element[Particles[P2]->element];
				double pplane = e -> pplane; 
				//cout<< "contact distance"<<Particles[P1]->h + pplane - dot (Particles[P2]->normal,	Particles[P1]->x)<<endl;
      				
				deltat_cont = ( Particles[P1]->h + pplane - dot (Particles[P2]->normal,	Particles[P1]->x) ) / (-delta_);								//Eq 3-142 
				Vec3_t Ri = Particles[P1]->x + deltat_cont * vr;	//Eq 3-139 Ray from SPH particle in the rel velocity direction

				//Check for contact in this time step 
				//Calculate time step for external forces
				double dt_fext = contact_force_factor * (Particles[P1]->Mass * 2. * norm(Particles[P1]->v) / norm (Particles[P1] -> contforce));	//Fraser 3-145
				// omp_set_lock(&Particles[P1]->my_lock);
        // Particles[P1] -> contforce = 0.; //RESET
				// omp_unset_lock(&Particles[P1]->my_lock);
				// if (dt_fext > deltat)
					// cout << "Time step size ("<<deltat<<" is larger than max allowable contact forces step ("<<dt_fext<<")"<<endl;
				if (deltat_cont < deltat){ //Originaly //	if (deltat_cont < std::min(deltat,dt_fext) 
					inside_time++;
					//cout << "Inside dt contact" <<endl;
					//Find point of contact Qj
					Vec3_t Qj = Particles[P1]->x + (Particles[P1]->v * deltat_cont) - ( Particles[P1]->h * Particles[P2]->normal); //Fraser 3-146
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
			
		// cout << "Particle Normal: "<<Particles[P2]->normal<<endl;
						// cout << "/////////////////////////////////////////" <<endl;
						// cout << " vr: "<< vr<<endl;
						//cout << "delta_: "<<delta_<<endl;
						// cout << "Particles[P1]->x"<< Particles[P1]->x<<endl;
						 // cout << "Particles[P2]->x"<< Particles[P2]->x<<endl;
						 // cout << "Particles[P2]->n"<< Particles[P2]->normal<<endl;
						 // cout << "Particles[P1]->h"<< Particles[P1]->h<<endl;
						 // cout << "pplane" << pplane<<endl;
						 // cout << "dot (Particles[P2]->normal,	Particles[P1]->x)" <<dot (Particles[P2]->normal,	Particles[P1]->x)<<endl;
						 //cout << "dt contact: "<<deltat_cont<<endl;
				
						//Recalculate vr (for large FEM mesh densities)
						//cout << "particle "<<P1 <<" inside element"<<endl;
						
						//Calculate penetration depth (Fraser 3-49)
						double delta = (deltat - deltat_cont) * delta_;
						// DEBUG THINGS, REMOVE ////////////
						inside_part[P1] ++;
						inside_part_count++;
						if (delta > max_delta) max_delta = delta;
						if (delta < min_delta) min_delta = delta;
					///////////////////////////
				
						//cout << "delta: "<<delta<<endl;
						// omp_set_lock(&Particles[P1]->my_lock);
						// if (Particles[P1]->delta_pl_strain > 0.) {
							// //cout << "recalc sitffness"<<endl;
							// double dS = pow(Particles[P1]->Mass/Particles[P1]->Density,0.33333); //Fraser 3-119;
							// //Particles[P1]-> cont_stiff = Particles[P1]->mat->Elastic().E()*0.01 * dS;
							// Particles[P1]-> cont_stiff = Particles [P1]-> Et_m * dS;
							// //cout << "recalculated: "<< Particles[P1]-> cont_stiff<<endl;
						// }
						// omp_unset_lock(&Particles[P1]->my_lock);
						
						// DAMPING
						//Calculate SPH and FEM elements stiffness (series)
						//Since FEM is assumed as rigid, stiffness is simply the SPH one 
						double kij = PFAC * Particles[P1]-> cont_stiff;
						double omega = sqrt (kij/Particles[P1]->Mass);
						double psi_cont = 2. * Particles[P1]->Mass * omega * DFAC; // Fraser Eqn 3-158

						omp_set_lock(&Particles[P1]->my_lock);
						Particles[P1] -> contforce = (kij * delta - psi_cont * delta_) * Particles[P2]->normal; // NORMAL DIRECTION, Fraser 3-159
            Particles[P1] -> delta_cont = delta;
						omp_unset_lock(&Particles[P1]->my_lock);
						force2 = dot(Particles[P1] -> contforce,Particles[P1] -> contforce);
						
						// TANGENTIAL COMPONENNT
						// Fraser Eqn 3-167
						// TODO - recalculate vr here too!
						Vec3_t tgvr, tgdir;
						double norm_tgvr;
						if (friction > 0.) {						
							if ( norm (vr)  != 0.0 ) {
								//TODO: THIS VELOCITY SHOULD BE THE CORRECTED ONE 
								//Vec3_t tgvr  = vr - dot(vr,Particles[P2]->normal) * Particles[P2]->normal;
								// Is Fraser thesis is explained better 
								Vec3_t tgvr = vr + delta_ * Particles[P2]->normal;  // -dot(vr,normal) * normal
								norm_tgvr = norm(tgvr);
								tgdir = tgvr / norm_tgvr;
								omp_set_lock(&Particles[P1]->my_lock);
								Particles[P1] -> tgdir = tgdir; // NORMAL DIRECTION, Fraser 3-159
								omp_unset_lock(&Particles[P1]->my_lock);
							}
						}
						
						// if (force2 > (1.e10))
							// Particles[P1] -> contforce = 1.e5;
						dt_fext = contact_force_factor * (Particles[P1]->Mass * 2. * norm(Particles[P1]->v) / norm (Particles[P1] -> contforce));

						if (dt_fext < min_force_ts_){
							min_force_ts_ = dt_fext;
							if (dt_fext > 0)
								this -> min_force_ts = min_force_ts_;
						}
						omp_set_lock(&Particles[P1]->my_lock);
						Particles[P1] -> a += Particles[P1] -> contforce / Particles[P1] -> Mass; 
						omp_unset_lock(&Particles[P1]->my_lock);
						//cout << "contforce "<<Particles[P1] -> contforce<<endl;
						
						Vec3_t tgforce;
						if (friction > 0.) {
							if (norm_tgvr > VMIN_FOR_FRICTION){

							// //TG DIRECTION
								tgforce = friction * norm(Particles[P1] -> contforce) * tgdir;
								omp_set_lock(&Particles[P1]->my_lock);
								Particles[P1] -> a -= tgforce / Particles[P1] -> Mass; 
								omp_unset_lock(&Particles[P1]->my_lock);
								//cout << "tg force "<< tgforce <<endl;
								
								if (cont_heat_gen) {
								//Fraser Eqns 3.109 to 3.111
								//lambda is fraction passed to 
								omp_set_lock(&Particles[P1]->my_lock);
								Particles[P1]->q_fric_work = dot(tgforce,vr); //J/(m3.s)
								omp_unset_lock(&Particles[P1]->my_lock);
								}
							
							}
						}
						
						if   (force2 > max_contact_force ) max_contact_force = force2;
						else if (force2 < min_contact_force ) min_contact_force = force2;
						inside_pairs++;
					}// if inside
				} //deltat <min

			}//delta_ > 0 : PARTICLES ARE APPROACHING EACH OTHER
		}//Contact Pairs
	}//Nproc
  //cout << "END CONTACT----------------------"<<endl;
	max_contact_force = sqrt (max_contact_force);
	//min_contact_force = sqrt (min_contact_force);
	//cout << "Inside pairs count: "<<inside_part_count<<", Inside time: "<<inside_time<<", Total cont pairs" << cont_pairs <<endl;
	inside_part_count = 0;
	int cont_force_count = 0;
	// for (int i = 0;i<Particles.Size();i++){
		// //DO THIS IN SERIAL (NOT PARALLEL) MODE OR BLOCK THIS IN PRAGMA
		// if (inside_part[i]>0)
			// inside_part_count++;
		// if (inside_part[i]>1)
			// cout << "WARNING, SAME PARTICLE HAS 2 DIFFERENT RIGID CONTACT POINTS. "<<endl;
		// // if (norm(Particles[i]->contforce)>0.)
			// // cont_force_count++;
	// }
	// cout << "Inside particle count: "<<inside_part_count<<endl;
	// cout << "Max penetration: "<<max_delta<<", min penetration: "<<min_delta<<endl; 
	//cout << "Particles with contact force: "<<cont_force_count<<endl;
	
	if (max_contact_force > 0.){
		//cout << "Min Contact Force"<< min_contact_force<<"Max Contact Force: "<< max_contact_force << "Time: " << Time << ", Pairs"<<inside_pairs<<endl;
		//cout << " Min tstep size: " << min_force_ts << ", current time step: " << deltat <<endl;
		//TEMP
		// if (min_force_ts> 0)
			// deltat = min_force_ts;
	}
	//Correct time step!
//	std::min(deltat,dt_fext)
}

}; //SPH