#include "matvec.h" 
#define VMIN_FOR_FRICTION	1.e-3

#define VMAX_FOR_STA_FRICTION	1.e-3

namespace SPH {

	
void Domain::AddTrimeshParticles(TriMesh *mesh, const float &hfac, const int &id){
	//if (meshcount==0)
    first_fem_particle_idx.push_back(Particles.Size());
	// else
    // first_fem_particle_idx.push_back(trimesh[meshcount-1]->element.Size());
  
	double Density =0.;
	double h;
	bool Fixed = false;	//Always are fixed ...
	contact_surf_id.push_back(id);
	//trimesh[m] = &mesh;
  trimesh.push_back(mesh);
	cout << "Adding particles"<<endl;
	for ( int e = 0; e < mesh->element.Size(); e++ ){
		Vec3_t pos = mesh->element[e]->centroid;
		h = hfac * mesh->element[e]->radius;
		Particles.Push(new Particle(id,pos,Vec3_t(0,0,0),0.0,Density,h,Fixed));
		Particles[first_fem_particle_idx[meshcount] + e] -> normal  = mesh->element[e] -> normal;
		Particles[first_fem_particle_idx[meshcount] + e] -> element = e; 
    Particles[first_fem_particle_idx[meshcount] + e] -> mesh = meshcount;
	}
	cout << Particles.Size() - first_fem_particle_idx[meshcount] << "particles added with ID " << contact_surf_id[meshcount] <<endl;
	cout << first_fem_particle_idx[meshcount] << " is the first solid particle index."<<endl;
  meshcount++;
}

//PARTICLES POSITIONS IS USED IN MOVE!
inline void Domain::UpdateContactParticles(){
  for (int m=0; m<trimesh.size();m++){
    for ( int e = 0; e < trimesh[m]->element.Size(); e++ ){
      Vec3_t v = 0.;
      for (int en = 0;en<3;en++)
        v += *trimesh[m] -> node_v[trimesh[m]->element[e] ->node[en]];
      Particles[first_fem_particle_idx[m] + e] -> v = 
      Particles[first_fem_particle_idx[m] + e] -> va = 
      Particles[first_fem_particle_idx[m] + e] -> vb = v/3.;
      Particles[first_fem_particle_idx[m] + e] -> a = 0.; 
      Particles[first_fem_particle_idx[m] + e] -> normal  = trimesh[m]->element[e] -> normal;
      //cout << "v "<< v/3.<<", n "<<Particles[first_fem_particle_idx + e] -> normal<<endl;
    } 
  }
}

inline void Domain::ContactNbSearch(){
	//cout << "Performing Nb Search"<<endl;
  		size_t P1,P2;
	cont_pairs = 0;
	#pragma omp parallel for schedule (static) private (P1,P2) num_threads(Nproc)
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
		//cout << "Rig Pair size"<<RIGPairs[k].Size()<<endl;
		for (size_t a=0; a<RIGPairs[k].Size();a++) {
			
			P1	= RIGPairs[k][a].first;
			P2	= RIGPairs[k][a].second;	
			// cout << "Fist Particle index, ID: "<<P1 << ", " << Particles[P1]->ID  <<endl;
			// cout << "Sec Particle index, ID: "<<P2 << ", " << Particles[P2]->ID  <<endl;
      bool is_contact = false;
      for (int m=0;m<meshcount;m++){
        if (Particles[P1]->ID == contact_surf_id[m] || Particles[P2]->ID == contact_surf_id[m] )
          is_contact = true;
      }
      if (is_contact) {
        if (Particles[P1]->ID == id_free_surf || Particles[P2]->ID == id_free_surf ) {
          Vec3_t xij	= Particles[P1]->x - Particles[P2]->x;

          //cout << "r, rcutoff, h1, h2"<< r << ", "<< rcutoff << ", "<< Particles[P1]->h <<", "<<Particles[P2]->h<<endl;
          if ( norm (Particles[P1]->x - Particles[P2]->x) < ( Particles[P1]->h + Particles[P2]->h ) ){ //2* cutoff being cutoff (h1+h2)/2
          //cout << "Found contact pair: "<< P1 << ", " << P2 << endl;
          //ContPairs[k].Push(std::make_pair(P1, P2));
          ContPairs[k].Push(RIGPairs[k][a]);
          cont_pairs++;
          //If the problem is not thermal (only mechanic)
          //Could be deleted this pair in Whole Pairs
          }	//belongs to free surf
        }
        // Pair is removed either way is inside cutoff radius or is or in!
        
        //RIGPairs[k].DelItem(a);//Erase, NOT EFFICIENT
        //a--;
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
  
  	#pragma omp parallel for schedule (dynamic) num_threads(Nproc)
    for (int i=0 ; i<Nproc ; i++)  //In the original version this was calculated after
		RIGPairs[i].Clear();
    
	// cout << "Contact Pairs Count"<<endl;
	// for (int k=0; k<Nproc;k++) 
		// cout << ContPairs[k].Size()<<", ";
	 // cout <<endl;
}

inline void Domain::CalcContactInitialGap(){
  cout << "Calculaint initial gap"<<endl;
	double min_delta,max_delta;
	min_delta = 1000.; max_delta = 0.;
	int inside_time;
	

	max_contact_force = 0.;
	double min_contact_force = 1000.;
	int inside_pairs = 0;
	double force2 = 0.;
  double distance;
  double crit;
  
  double kij, omega,psi_cont;
  int i,j;  //For inside testing
	
	int P1,P2;
  int m;

  //Vec3_t vr[Particles.Size()];
  Vec3_t vr;

  Element* e;

  Vec3_t atg;
  double mindist = 1000.;
  double maxdist = -1000.;
  double delta_;
  
	#pragma omp parallel for schedule (static) private(P1,P2,vr,delta_,m,distance, e) num_threads(Nproc)
  //tgforce
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
		//cout << "Contact pair size: "<<RIGPairs[k].Size()<<endl;
		for (size_t a = 0; a < RIGPairs[k].Size();a++) {
			//P1 is SPH particle, P2 is CONTACT SURFACE (FEM) Particle
      for (int m=0;m<meshcount;m++){
        if (Particles[RIGPairs[k][a].first]->ID == contact_surf_id[m] ) 	{ 	//Cont Surf is partcicles from FEM
          P1 = RIGPairs[k][a].second; P2 = RIGPairs[k][a].first; 	}
        else {
          P1 = RIGPairs[k][a].first; P2 = RIGPairs[k][a].second; } 
      } 
			vr = Particles[P1]->v - Particles[P2]->v;		//Fraser 3-137
			delta_ = - dot( Particles[P2]->normal , vr);	//Penetration rate, Fraser 3-138
      
      //cout << "p1 vel "<<Particles[P1]->v << "p2 vel "<< Particles[P2]->v <<endl;
      // cout << "distance "<< Particles[P1]->x - Particles[P2]->x<<endl;
			//Check if SPH and fem particles are approaching each other
			if (delta_ > 0 ){
        m = Particles[P2]->mesh;
        //cout << "particle Mesh "<< m<<endl;
        //if (m!=-1){
        e = trimesh[m]-> element[Particles[P2]->element];
              
        distance = -( Particles[P1]->h + trimesh[m]-> element[Particles[P2]->element] -> pplane 
                      - dot (Particles[P2]->normal,	Particles[P1]->x) ) ;								//Eq 3-142 
        //cout << "pplane: "<<trimesh-> element[Particles[P2]->element] -> pplane <<endl;        
        if (distance  < mindist){
          omp_set_lock(&dom_lock);
            mindist = distance;
          omp_unset_lock(&dom_lock);
        } else if (distance  > maxdist){
          omp_set_lock(&dom_lock);
            maxdist = distance;
          omp_unset_lock(&dom_lock);        
        }
      }
    }
  }
    //cout << "Min contact gap is " << mindist<<endl;
    //cout << "Max contact gap is " << maxdist<<endl;    

}

//////////////////////////////// 
//// 
////////////////////////////////
inline void Domain::CalcContactForces(){
	
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
	int inside_time,inside_geom;
	
  //#pragma omp parallel for num_threads(Nproc)
  for (int i = 0;i<Particles.Size();i++){
		omp_set_lock(&Particles[i]->my_lock);
		Particles[i] -> contforce = 0.; //RESET    
		Particles[i] -> delta_cont = 0.; //RESET    
		Particles[i] -> tgdir = 0.;				//TODO: DELETE (DEBUG) 
		omp_unset_lock(&Particles[i]->my_lock);
		inside_part[i] = 0;
		inside_time=inside_geom=0;
  }
 
	
	max_contact_force = 0.;
	double min_contact_force = 1000.;
	int inside_pairs = 0;
	double force2 = 0.;
	double delta_ = 0.;
	double deltat_cont;
  double crit;
  
  double kij, omega,psi_cont;
  int i,j;  //For inside testing
	
	int P1,P2;
  Vec3_t tgforce;
  Vec3_t Qj[Particles.Size()]; //Things not allowed
  //Vec3_t vr[Particles.Size()];
  Vec3_t vr;
  double dt_fext;
  Element* e;
  bool inside;
  
  Vec3_t tgvr, tgdir;
  double norm_tgvr;
  double max_vr = 0.;
  int m;
 
  Vec3_t atg;
  
  contact_force_sum = 0.;
  
  int max_reached_part = 0; //TEST
  int sta_frict_particles = 0;
	#pragma omp parallel for schedule (static) private(P1,P2,vr,delta_,deltat_cont, m, inside,i,j,crit,force2,dt_fext,kij,omega,psi_cont,e,tgforce,tgvr,norm_tgvr,tgdir,atg) num_threads(Nproc)
  //tgforce
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
      bool is_first = false;
      for (int m=0;m<meshcount;m++){
        if (Particles[ContPairs[k][a].first]->ID == contact_surf_id[m])
          is_first = true;
      }
      

      if (is_first) 	{ 	//Cont Surf is partcicles from FEM
        P1 = ContPairs[k][a].second; P2 = ContPairs[k][a].first; 	}
      else {
        P1 = ContPairs[k][a].first; P2 = ContPairs[k][a].second; } 
      
			vr = Particles[P1]->v - Particles[P2]->v;		//Fraser 3-137
			//cout << "Particle P1v: "<<Particles[P1]->v<<endl;
			//cout << "Particle P2v: "<<Particles[P2]->v<<endl;
			//ok, FEM Particles normals can be calculated by two ways, the one used to
			//calculate SPH ones, and could be given by mesh input
			//delta_ Is the projection of relative velocity 
			delta_ = - dot( Particles[P2]->normal , vr);	//Penetration rate, Fraser 3-138
      
      //cout << "p1 vel "<<Particles[P1]->v << "p2 vel "<< Particles[P2]->v <<endl;
      // cout << "distance "<< Particles[P1]->x - Particles[P2]->x<<endl;
			//Check if SPH and fem particles are approaching each other
			if (delta_ > 0 ){
        m = Particles[P2]->mesh;
        //cout << "particle Mesh "<< m<<", " <<"particle " << P2<<endl;
				e = trimesh[m]-> element[Particles[P2]->element];
				//double pplane = trimesh-> element[Particles[P2]->element] -> pplane; 
				//cout<< "contact distance"<<Particles[P1]->h + trimesh-> element[Particles[P2]->element] -> pplane - dot (Particles[P2]->normal,	Particles[P1]->x)<<endl;
      				
				deltat_cont = ( Particles[P1]->h + trimesh[m]-> element[Particles[P2]->element] -> pplane 
                      - dot (Particles[P2]->normal,	Particles[P1]->x) ) / (- delta_);								//Eq 3-142 
				//Vec3_t Ri = Particles[P1]->x + deltat_cont * vr;	//Eq 3-139 Ray from SPH particle in the rel velocity direction

				//Check for contact in this time step 
				//Calculate time step for external forces
				// dt_fext = contact_force_factor * (Particles[P1]->Mass * 2. * norm(Particles[P1]->v) / norm (Particles[P1] -> contforce));	//Fraser 3-145
				// omp_set_lock(&Particles[P1]->my_lock);
        // Particles[P1] -> contforce = 0.; //RESET
				// omp_unset_lock(&Particles[P1]->my_lock);
				// if (dt_fext > deltat)
					// cout << "Time step size ("<<deltat<<" is larger than max allowable contact forces step ("<<dt_fext<<")"<<endl;
        //cout << "deltat_cont "<<deltat_cont<< "deltat" << deltat <<endl; 
				if (deltat_cont < deltat){ //Originaly //	if (deltat_cont < std::min(deltat,dt_fext) 
					inside_time++;
					//cout << "Inside dt contact" <<endl;
					//Find point of contact Qj
					Qj[P1] = Particles[P1]->x + (Particles[P1]->v * deltat_cont) - ( Particles[P1]->h * Particles[P2]->normal); //Fraser 3-146
					//Check if it is inside triangular element
					//Find a vector 
					//Fraser 3-147
					inside = true;
					i=0;		
					while (i<3 && inside){
						j = i+1;	if (j>2) j = 0;
						crit = dot (cross ( *trimesh[m]->node[e -> node[j]] 
                                        - *trimesh[m]->node[e -> node[i]],
                                        Qj[P1]  - *trimesh[m]->node[e -> node[i]]),
															Particles[P2]->normal);
						if (crit < 0.0) inside = false;
						i++;
					}
					
					if (inside ) { //Contact point inside element, contact proceeds
            inside_geom++;
            //cout << "Particle Normal: "<<Particles[P2]->normal<<endl;
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
						// if (delta > max_delta) max_delta = delta;
						// if (delta < min_delta) min_delta = delta;
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
						kij = PFAC * Particles[P1]-> cont_stiff;
						omega = sqrt (kij/Particles[P1]->Mass);
						psi_cont = 2. * Particles[P1]->Mass * omega * DFAC; // Fraser Eqn 3-158

						omp_set_lock(&Particles[P1]->my_lock);
						Particles[P1] -> contforce = (kij * delta - psi_cont * delta_) * Particles[P2]->normal; // NORMAL DIRECTION, Fraser 3-159
            Particles[P1] -> delta_cont = delta;
						omp_unset_lock(&Particles[P1]->my_lock);

            omp_set_lock(&dom_lock);            
              contact_force_sum += norm(Particles[P1] ->contforce);
            omp_unset_lock(&dom_lock);	
            
						force2 = dot(Particles[P1] -> contforce,Particles[P1] -> contforce);
						
						// TANGENTIAL COMPONENNT DIRECTION
						// Fraser Eqn 3-167
						// TODO - recalculate vr here too!
            
            //removed if, this is calculated always
            tgvr = vr + delta_ * Particles[P2]->normal;  // -dot(vr,normal) * normal, FRASER 3-168
            norm_tgvr = norm(tgvr);  
            tgdir = tgvr / norm_tgvr;              

            if (friction_sta > 0.) { 
              atg = Particles[P1] -> a - dot (Particles[P1] -> a,Particles[P2]->normal)*Particles[P2]->normal;
              if ( (norm(atg) * Particles[P1] -> Mass) < friction_sta * norm(Particles[P1] -> contforce)
                && norm_tgvr < VMAX_FOR_STA_FRICTION) {
                omp_set_lock(&Particles[P1]->my_lock);
                //cout << "particle accel x and y before "<<Particles[P1] -> a[0]<<", "<<Particles[P1] -> a[1] <<endl;
                //cout << "atg "<<atg<<endl;
                Particles[P1] -> tgdir = atg;
                Particles[P1] -> a -= atg; 
                
                // THIS CRASHES
                Particles[P1] -> v = Particles[P1] -> va = Particles[P1] -> vb = Particles[P2] -> v; 
                
                
                //cout << "particle 2 vel "<<Particles[P2] -> v<<endl;
                //cout << "particle accel x and y after"<<Particles[P1] -> a[0]<<", "<<Particles[P1] -> a[1] <<endl;
                //cout << "particle vx vy "<< Particles[P1] -> v[0]<<", "<<Particles[P1] -> a[1] <<endl;
                omp_unset_lock(&Particles[P1]->my_lock);
                // omp_set_lock(&dom_lock);
                    // sta_frict_particles++;
                // omp_unset_lock(&dom_lock);
              }
              else {
                // omp_set_lock(&dom_lock);
                  // max_reached_part++;
                // omp_unset_lock(&dom_lock);
                //cout << "Max force reached! particle vx vy "<< Particles[P1] -> v[0]<<", "<<Particles[P1] -> v[1] <<endl;
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
						
						if (friction_dyn > 0.) {
              //if (fric_type==Fr_Dyn){
                if (norm_tgvr > /*0.0*/VMIN_FOR_FRICTION){

                // //TG DIRECTION
                  tgforce = friction_dyn * norm(Particles[P1] -> contforce) * tgdir;
                  omp_set_lock(&Particles[P1]->my_lock);
                  Particles[P1] -> a -= tgforce / Particles[P1] -> Mass; 
                  omp_unset_lock(&Particles[P1]->my_lock);
                  //cout << "tg force "<< tgforce <<endl;
                  
                  if (cont_heat_gen) {
                  //Fraser Eqns 3.109 to 3.111
                  //lambda is fraction passed to 
                  omp_set_lock(&Particles[P1]->my_lock);
                  Particles[P1]->q_fric_work = dot(tgforce,vr); //J/(m3.s)
                  //cout<< Particles[P1]->q_fric_work<<endl;
                  omp_unset_lock(&Particles[P1]->my_lock);
                  }
                
                }
              //}
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
    //cout << "particles surpassed max fr force"<<max_reached_part<< ", below force: " <<sta_frict_particles<<endl;
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