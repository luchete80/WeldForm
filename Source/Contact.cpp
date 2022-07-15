#include "matvec.h" 
#define VMIN_FOR_FRICTION	1.e-3

#define VMAX_FOR_STA_FRICTION	1.e-3
#include "Plane.h"

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
  mesh->id = id;
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
  double normal_cf;
  Vec3_t atg;
  bool end;
  contact_force_sum = 0.;
  double delta;
  Vec3_t delta_tg;
  int max_reached_part = 0; //TEST
  int sta_frict_particles = 0;
  int stra_restr = 0; //restricted static
	#pragma omp parallel for schedule (static) private(P1,P2,end,vr,delta_,delta, deltat_cont, m, inside,delta_tg,i,j,crit,normal_cf,force2,dt_fext,kij,omega,psi_cont,e,tgforce,tgvr,norm_tgvr,tgdir,atg) num_threads(Nproc)
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
    // TODO: DO THIS ONCE PER PARTICLE
    // int a = 0;
    // end = false;
    // while (!end){ //Per particle
    
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
      	//  Distance plane to point: 
        //xi . n - d = pplane  (Point at side of normal)	
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
            end=true;
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
						delta = (deltat - deltat_cont) * delta_;
						// DEBUG THINGS, REMOVE ////////////
						inside_part[P1] ++;
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
            
            //psi_cont = Particles[i]->Cs *Particles[i]->Density;
            //cout << "psi_cont "<<psi_cont/DFAC<<endl; 
            normal_cf = (kij * delta - psi_cont * delta_);
						omp_set_lock(&Particles[P1]->my_lock);
						Particles[P1] -> contforce = normal_cf * Particles[P2]->normal; // NORMAL DIRECTION, Fraser 3-159
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
            atg = Particles[P1] -> a - dot (Particles[P1] -> a,Particles[P2]->normal)*Particles[P2]->normal;
            //ONCE SAVED atg, now can change acceleration by contact!

						// if (force2 > (1.e10))
							// Particles[P1] -> contforce = 1.e5;
						dt_fext = contact_force_factor * (Particles[P1]->Mass * 2. * norm(Particles[P1]->v) / norm (Particles[P1] -> contforce));

						if (dt_fext < min_force_ts_){
							min_force_ts_ = dt_fext;
							if (dt_fext > 0)
								this -> min_force_ts = min_force_ts_;
						}
            //Acceleration is summed up on integration
						// omp_set_lock(&Particles[P1]->my_lock);
						// Particles[P1] -> a += Particles[P1] -> contforce / Particles[P1] -> Mass; 
						// omp_unset_lock(&Particles[P1]->my_lock);
						//cout << "contforce "<<Particles[P1] -> contforce<<endl;

            
            if (fric_type == Fr_Bound){
                omp_set_lock(&Particles[P1]->my_lock);
                Particles[P1] -> a -= atg; 
                //cout << "after adj: accel: "<< Particles[P1] -> a<<endl;
                //cout << "atg: "<< atg<<endl;
                omp_unset_lock(&Particles[P1]->my_lock);
                
              
            } else {
              if (friction_sta > 0.) { 
                //delta_tg = -vr * (deltat - deltat_cont) - ( delta * Particles[P2]->normal);  //THIS IS OPPOSITE TO DIRECTION
                
                if (P1 == 12415){
                  //CONTROL, particle 12415x -0.0075, y 0.1275, z 0.604
                cout << "delta tg 1 "<<delta_tg<<endl;
                delta_tg = (vr - dot(vr,Particles[P2]->normal)*Particles[P2]->normal)* (deltat - deltat_cont); //Viewed from P1
                tgforce = (kij * delta_tg - psi_cont * delta_);
                cout << "delta tg 2 "<<delta_tg<<endl;
                }
                
                if (norm(tgforce) < friction_sta * normal_cf ){
                  omp_set_lock(&Particles[P1]->my_lock);
                    Particles[P1] -> contforce += tgforce;
                    Particles[P1] -> a -= tgforce / Particles[P1]->Mass;  // //Eqn 30. Zhan
                  omp_unset_lock(&Particles[P1]->my_lock);
                }
                // if ( (norm(atg) * Particles[P1] -> Mass) < friction_sta * norm(Particles[P1] -> contforce)
                  // && norm_tgvr < VMAX_FOR_STA_FRICTION) {
                  // omp_set_lock(&Particles[P1]->my_lock);
                  // //cout << "particle accel x and y before "<<Particles[P1] -> a[0]<<", "<<Particles[P1] -> a[1] <<endl;
                  // //cout << "atg "<<atg<<endl;
                  // Particles[P1] -> tgdir = atg;
                  // Particles[P1] -> a -= atg; 
                  // stra_restr++;
                  // // THIS CRASHES
                  // omp_unset_lock(&Particles[P1]->my_lock);
                  // // omp_set_lock(&dom_lock);
                      // // sta_frict_particles++;
                  // // omp_unset_lock(&dom_lock);
                  
                // }
                
              }
              
              // if (friction_dyn > 0.) {
                // //if (fric_type==Fr_Dyn){
                  // if (norm_tgvr > /*0.0*/VMIN_FOR_FRICTION){

                  // // //TG DIRECTION
                    // tgforce = friction_dyn * norm(Particles[P1] -> contforce) * tgdir;
                    // omp_set_lock(&Particles[P1]->my_lock);
                    // Particles[P1] -> a -= tgforce / Particles[P1] -> Mass; 
                    // omp_unset_lock(&Particles[P1]->my_lock);
                    // //cout << "tg force "<< tgforce <<endl;
                    
                    // if (cont_heat_gen) {
                    // //Fraser Eqns 3.109 to 3.111
                    // //lambda is fraction passed to 
                    // omp_set_lock(&Particles[P1]->my_lock);
                    // Particles[P1]->q_fric_work = dot(tgforce,vr); //J/(m3.s)
                    // //cout<< Particles[P1]->q_fric_work<<endl;
                    // omp_unset_lock(&Particles[P1]->my_lock);
                    // }
                  
                  // }
                // //}
              // }
						}//Friction type
						if   (force2 > max_contact_force ) max_contact_force = force2;
						else if (force2 < min_contact_force ) min_contact_force = force2;
						inside_pairs++;
            
					}// if inside
				} //deltat <min

			}//delta_ > 0 : PARTICLES ARE APPROACHING EACH OTHER
		
      // a++;
      // if (a==ContPairs[k].Size())
        // end=true;
    }//Contact Pairs
	}//Nproc
  //cout << "END CONTACT----------------------"<<endl;
	max_contact_force = sqrt (max_contact_force);
  //cout << "contact_force_sum "<<contact_force_sum<<endl;
	//min_contact_force = sqrt (min_contact_force);
  //cout << "Inside pairs count: "<<inside_geom<<", Inside time: "<<inside_time<<", statically restricted " << stra_restr<<endl;
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
	
	//if (max_contact_force > 0.){
    //cout << "particles surpassed max fr force"<<max_reached_part<< ", below force: " <<sta_frict_particles<<endl;
		//cout << "Min Contact Force"<< min_contact_force<<"Max Contact Force: "<< max_contact_force << "Time: " << Time << ", Pairs"<<inside_pairs<<endl;
		//cout << " Min tstep size: " << min_force_ts << ", current time step: " << deltat <<endl;
		//TEMP
		// if (min_force_ts> 0)
			// deltat = min_force_ts;
	//}
	//Correct time step!
//	std::min(deltat,dt_fext)
}


//THIS IS THE SEO2006 METHOD, WHICH CHECKS ONLY PENETRATION (NON VELOCITY DIRECTION)
inline void Domain::CalcContactForces2(){

	double min_force_ts_=1000.;
// https://stackoverflow.com/questions/10850155/whats-the-difference-between-static-and-dynamic-schedule-in-openmp

	////////////////////////
	// DEBUG THINGS //////
	int inside_part[Particles.Size()];
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
	double delta_ = 0.;
  double delta;
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
  Vec3_t delta_tg;
  double max_vr = 0.;
  int m;
  double normal_cf;
  Vec3_t du;// If no contact
  
  bool ref_accel = true;
  
  Vec3_t ref_tg;
 
  Vec3_t atg;
  bool end;
  contact_force_sum = 0.;
  double dist;
  int max_reached_part = 0; //TEST
  int sta_frict_particles = 0;
  int stra_restr = 0; //restricted static
	Vec3_t x_pred;
  #pragma omp parallel for schedule (static) private(P1,P2,end,vr,dist, delta_tg, delta_,delta, du, normal_cf, ref_tg, m, x_pred, inside,i,j,crit,dt_fext,kij,omega,psi_cont,e,tgforce,tgvr,norm_tgvr,tgdir,atg) num_threads(Nproc)
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
    // TODO: DO THIS ONCE PER PARTICLE
    // int a = 0;
    // end = false;
    // while (!end){ //Per particle
    
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
			//if (delta_ > 0 ){
        m = Particles[P2]->mesh;
        //cout << "particle Mesh "<< m<<", " <<"particle " << P2<<endl;
				e = trimesh[m]-> element[Particles[P2]->element];

        dist =  dot (Particles[P2]->normal, Particles[P1]->x ) - trimesh[m]-> element[Particles[P2]->element] -> pplane;
        if( dist  < Particles[P1]->h) {

          Qj[P1] = Particles[P1]->x - dist * Particles[P2]->normal;
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
            end=true;

            delta = Particles[P1]->h - dist;
            //cout << "dist "<<dist<<", h "<<Particles[P1]->h<< ", delta "<<delta<<endl;
						// DAMPING
						//Calculate SPH and FEM elements stiffness (series)
						//Since FEM is assumed as rigid, stiffness is simply the SPH one 
						kij = PFAC * Particles[P1]-> cont_stiff;
						omega = sqrt (kij/Particles[P1]->Mass);
						psi_cont = 2. * Particles[P1]->Mass * omega * DFAC; // Fraser Eqn 3-158 
            //normal_cf = 2.0 * Particles[P1]->Mass /(deltat*deltat )*delta;

						omp_set_lock(&Particles[P1]->my_lock);
						//Particles[P1] -> contforce = normal_cf *  Particles[P2]->normal; 
            Particles[P1] -> contforce = (kij * delta - psi_cont * delta_) * Particles[P2]->normal; // NORMAL DIRECTION, Fraser 3-159
            Particles[P1] -> delta_cont = delta;
						omp_unset_lock(&Particles[P1]->my_lock);

            omp_set_lock(&dom_lock);            
              contact_force_sum += norm(Particles[P1] ->contforce);
            omp_unset_lock(&dom_lock);	
            //inside_pairs++;
						
						// TANGENTIAL COMPONENNT DIRECTION
						// Fraser Eqn 3-167
						// TODO - recalculate vr here too!
            
            //removed if, this is calculated always
            tgvr = vr + delta_ * Particles[P2]->normal;  // -dot(vr,normal) * normal, FRASER 3-168
            norm_tgvr = norm(tgvr);  
            tgdir = tgvr / norm_tgvr;              
            atg = Particles[P1] -> a - dot (Particles[P1] -> a,Particles[P2]->normal)*Particles[P2]->normal;
            //ONCE SAVED atg, now can change acceleration by contact!

						dt_fext = contact_force_factor * (Particles[P1]->Mass * 2. * norm(Particles[P1]->v) / norm (Particles[P1] -> contforce));

						if (dt_fext < min_force_ts_){
							min_force_ts_ = dt_fext;
							if (dt_fext > 0)
								this -> min_force_ts = min_force_ts_;
						}
						// omp_set_lock(&Particles[P1]->my_lock);
						// Particles[P1] -> a += Particles[P1] -> contforce / Particles[P1] -> Mass; 
						// omp_unset_lock(&Particles[P1]->my_lock);
						//cout << "contforce "<<Particles[P1] -> contforce<<endl;
            
            // TANGENTIAL BEHAVIOUR
            
            //Wang2013, but applied to current step
            //ORIGINAL VERSION, BASED IN CURRENT STEP
            du = vr * deltat;

            //Wang2013, but applied to current step
            // x_pred = Particles[P1]->x + Particles[P1]->v * deltat + Particles[P1]->a * deltat * deltat/2.0;
            // du = x_pred - Particles[P1] ->x - Particles[P2] -> v * deltat ;
            
            delta_tg = du - dot(du, Particles[P2]->normal)*Particles[P2]->normal;
            tgforce = kij * delta_tg;
            //tgforce = (kij * delta_tg - psi_cont * delta_);
                
           if (ref_accel) ref_tg = atg * Particles[P1]->Mass;
           else           ref_tg = tgforce;
               
            if (friction_sta > 0.) { 
                // //delta_tg = -vr * (deltat - deltat_cont) - ( delta * Particles[P2]->normal);  //THIS IS OPPOSITE TO DIRECTION
                
                if (P1 == 11312){
                  //CONTROL, particle 12415x -0.0075, y 0.1275, z 0.604
                cout << "delta tg 2 "<<delta_tg<<", delta "<<delta<<endl;
                cout << "normal du "<<dot(Particles[P1]->x_prev + vr, Particles[P2]->normal)*Particles[P2]->normal<<endl;
                cout << "tgforce " <<norm(tgforce) << ", mu N "<<friction_sta * norm(Particles[P1] -> contforce)<<endl;
                cout << "norm atg * mass : "<<norm(atg) * Particles[P1]->Mass<<endl;
                }
                
                if (norm(ref_tg) < friction_sta * norm(Particles[P1] -> contforce) ){
                  omp_set_lock(&Particles[P1]->my_lock);
                    Particles[P1] -> contforce -= tgforce / Particles[P1]->Mass; 
                    //Particles[P1] -> a -= tgforce / Particles[P1]->Mass;  // //Eqn 30. Zhan
                  omp_unset_lock(&Particles[P1]->my_lock);
                }
            }
            // if (fric_type == Fr_Bound){
                // omp_set_lock(&Particles[P1]->my_lock);
                // Particles[P1] -> a -= atg; 
                // Particles[P1] -> impose_vel = true;
                // Particles[P1] -> vc = Particles[P2]->vc; //Only tg vel?

                // omp_unset_lock(&Particles[P1]->my_lock);
                
              
            // } else {
              // //ORIGINAL
             // // if (friction > 0.) { 
                // // //if ( (norm(atg) * Particles[P1] -> Mass) < friction * norm(Particles[P1] -> contforce)) {
                  // // // cout << "acc " <<norm(atg) * Particles[P1] -> Mass<<endl;
                  // // // cout << "tg cont f "<<friction_sta * norm(Particles[P1] -> contforce)<<endl;
                  // // tgforce = friction * norm(Particles[P1] -> contforce) * tgdir;
                  // // omp_set_lock(&Particles[P1]->my_lock);
                    // // Particles[P1] -> tgdir = atg;
                    // // Particles[P1] -> a -= tgforce / Particles[P1] -> Mass; 
                    // // stra_restr++;
                  // // // cout << "particle accel x and y before "<<Particles[P1] -> a[0]<<", "<<Particles[P1] -> a[1] <<endl;
                  // // // cout << "atg "<<atg<<endl;
                  // // omp_unset_lock(&Particles[P1]->my_lock);
                  
                // // //}               
              // // }
              
             // if (friction_sta > 0.) { 
                // if ( (norm(atg) * Particles[P1] -> Mass) < friction_sta * norm(Particles[P1] -> contforce)
                  // && norm_tgvr < VMAX_FOR_STA_FRICTION) {
                  // // cout << "acc " <<norm(atg) * Particles[P1] -> Mass<<endl;
                  // // cout << "tg cont f "<<friction_sta * norm(Particles[P1] -> contforce)<<endl;
                  // omp_set_lock(&Particles[P1]->my_lock);
                    // Particles[P1] -> tgdir = atg;
                    // Particles[P1] -> a -= atg; 
                    // stra_restr++;
                  // // cout << "particle accel x and y before "<<Particles[P1] -> a[0]<<", "<<Particles[P1] -> a[1] <<endl;
                  // // cout << "atg "<<atg<<endl;
                  // omp_unset_lock(&Particles[P1]->my_lock);
                  
                // }               
              // }
              
              // if (friction_dyn > 0.) {
                // //if (fric_type==Fr_Dyn){
                  // if (norm_tgvr > /*0.0*/VMIN_FOR_FRICTION){

                  // // //TG DIRECTION
                    // tgforce = friction_dyn * norm(Particles[P1] -> contforce) * tgdir;
                    // omp_set_lock(&Particles[P1]->my_lock);
                    // Particles[P1] -> a -= tgforce / Particles[P1] -> Mass; 
                    // omp_unset_lock(&Particles[P1]->my_lock);
                    // //cout << "tg force "<< tgforce <<endl;
                    
                    // if (cont_heat_gen) {
                    // //Fraser Eqns 3.109 to 3.111
                    // //lambda is fraction passed to 
                    // omp_set_lock(&Particles[P1]->my_lock);
                    // Particles[P1]->q_fric_work = dot(tgforce,vr); //J/(m3.s)
                    // //cout<< Particles[P1]->q_fric_work<<endl;
                    // omp_unset_lock(&Particles[P1]->my_lock);
                    // }
                  
                  // }
                // //}
              // }              
              
            // }//!else Fr Bound
            
					}// if inside
        } //If distance is less than h
    }//Contact Pairs
	}//Nproc
  //cout << "END CONTACT----------------------"<<endl;
	max_contact_force = sqrt (max_contact_force);
	//min_contact_force = sqrt (min_contact_force);
	//cout << "Inside pairs count: "<<inside_geom<<", Inside time: "<<inside_time<<", statically restricted " << stra_restr<<endl;
	int cont_force_count = 0;
	
	//if (max_contact_force > 0.){
    //cout << "particles surpassed max fr force"<<max_reached_part<< ", below force: " <<sta_frict_particles<<endl;
		//cout << "Min Contact Force"<< min_contact_force<<"Max Contact Force: "<< max_contact_force << "Time: " << Time << ", Pairs"<<inside_pairs<<endl;
		//cout << " Min tstep size: " << min_force_ts << ", current time step: " << deltat <<endl;
		//TEMP
		// if (min_force_ts> 0)
			// deltat = min_force_ts;
	//}
	//Correct time step!
//	std::min(deltat,dt_fext)
}

// FOR ANALYTIC PLANES
inline void Domain::CalcContactForcesAnalytic(){

	double min_force_ts_=1000.;
// https://stackoverflow.com/questions/10850155/whats-the-difference-between-static-and-dynamic-schedule-in-openmp

	////////////////////////
	// DEBUG THINGS //////
	int inside_part[Particles.Size()];
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
  bool end;
  contact_force_sum = 0.;
  double delta;
  
  int max_reached_part = 0; //TEST
  int sta_frict_particles = 0;
  int stra_restr = 0; //restricted static
  #pragma omp parallel for schedule (static) num_threads(Nproc)
  for (int p=0; p<Particles.Size(); p++) {
		Vec3_t xij;
		double h,K;
		// Summing the smoothed pressure, velocity and stress for fixed particles from neighbour particles
		//IT IS CONVENIENT TO FIX SINCE FSMPairs are significantly smaller
		//cout << "Contact pair size: "<<ContPairs[k].Size()<<endl;

      
    vr = Particles[p]->v - Particles[P2]->v;		//Fraser 3-137
    delta_ = - dot( Particles[P2]->normal , vr);	//Penetration rate, Fraser 3-138
    
    //cout << "p vel "<<Particles[p]->v << "p2 vel "<< Particles[P2]->v <<endl;
    // cout << "distance "<< Particles[p]->x - Particles[P2]->x<<endl;
    //Check if SPH and fem particles are approaching each other
    if (delta_ > 0 ){
      m = Particles[P2]->mesh;
      //cout << "particle Mesh "<< m<<", " <<"particle " << P2<<endl;
      e = trimesh[m]-> element[Particles[P2]->element];
      //double pplane = trimesh-> element[Particles[P2]->element] -> pplane; 
      //cout<< "contact distance"<<Particles[p]->h + trimesh-> element[Particles[P2]->element] -> pplane - dot (Particles[P2]->normal,	Particles[p]->x)<<endl;
            
      deltat_cont = ( Particles[p]->h + trimesh[m]-> element[Particles[P2]->element] -> pplane 
                    - dot (Particles[P2]->normal,	Particles[p]->x) ) / (- delta_);								//Eq 3-142 

      if (deltat_cont < deltat){ //Originaly //	if (deltat_cont < std::min(deltat,dt_fext) 
        inside_time++;
        //cout << "Inside dt contact" <<endl;
        //Find point of contact Qj
        Qj[p] = Particles[p]->x + (Particles[p]->v * deltat_cont) - ( Particles[p]->h * Particles[P2]->normal); //Fraser 3-146
        
        //Calculate penetration depth (Fraser 3-49)
        delta = (deltat - deltat_cont) * delta_;

        // DAMPING
        //Calculate SPH and FEM elements stiffness (series)
        //Since FEM is assumed as rigid, stiffness is simply the SPH one 
        kij = PFAC * Particles[p]-> cont_stiff;
        omega = sqrt (kij/Particles[p]->Mass);
        psi_cont = 2. * Particles[p]->Mass * omega * DFAC; // Fraser Eqn 3-158

        omp_set_lock(&Particles[p]->my_lock);
        Particles[p] -> contforce = (kij * delta - psi_cont * delta_) * Particles[P2]->normal; // NORMAL DIRECTION, Fraser 3-159
        Particles[p] -> delta_cont = delta;
        omp_unset_lock(&Particles[p]->my_lock);
        
        contact_force_sum += norm(Particles[p] ->contforce);
        
        force2 = dot(Particles[p] -> contforce,Particles[p] -> contforce);
        
        // TANGENTIAL COMPONENNT DIRECTION
        // Fraser Eqn 3-167
        // TODO - recalculate vr here too!
        tgvr = vr + delta_ * Particles[P2]->normal;  // -dot(vr,normal) * normal, FRASER 3-168
        norm_tgvr = norm(tgvr);  
        tgdir = tgvr / norm_tgvr;              
        atg = Particles[p] -> a - dot (Particles[p] -> a,Particles[P2]->normal)*Particles[P2]->normal;
        //ONCE SAVED atg, now can change acceleration by contact!

        // if (force2 > (1.e10))
          // Particles[p] -> contforce = 1.e5;
        dt_fext = contact_force_factor * (Particles[p]->Mass * 2. * norm(Particles[p]->v) / norm (Particles[p] -> contforce));

        if (dt_fext < min_force_ts_){
          min_force_ts_ = dt_fext;
          if (dt_fext > 0)
            this -> min_force_ts = min_force_ts_;
        }
        Particles[p] -> a += Particles[p] -> contforce / Particles[p] -> Mass; 
        
        if (fric_type == Fr_Bound){
            Particles[p] -> a -= atg;             
        }//Friction type

      } //deltat <min

    }//delta_ > 0 : PARTICLES ARE APPROACHING EACH OTHER
  
    // a++;
    // if (a==ContPairs[k].Size())
      // end=true;

	}//for particles
  //cout << "END CONTACT----------------------"<<endl;
	max_contact_force = sqrt (max_contact_force);
  //cout << "contact_force_sum "<<contact_force_sum<<endl;
	//min_contact_force = sqrt (min_contact_force);
  //cout << "Inside pairs count: "<<inside_geom<<", Inside time: "<<inside_time<<", statically restricted " << stra_restr<<endl;
	int cont_force_count = 0;

//	std::min(deltat,dt_fext)
}


}; //SPH

#include "Contact_Wang.cpp"
