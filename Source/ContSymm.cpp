#include <vector>

namespace SPH {
  //////////////////////////////// 
//// From Zhan: A SPH framework for dynamic interaction between soil and
////            rigid body system with hybrid contact method
//// HERE PREDICTED VELOCITY IS PREV VELOCITY
////////////////////////////////
inline void Domain::CalcContactForcesWang(std::vector<int> * part_id, const Plane &plane){

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
  Vec3_t imp_force;
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
 
  //USE TRUE: FALSE DOES NOT WORK
  bool ref_accel = true;   //true: tg force is compared to current tg accel
  
  Vec3_t atg;
  bool end;
  contact_force_sum = 0.;
  double dist;
  int max_reached_part = 0; //TEST
  int sta_frict_particles = 0;
  int stra_restr = 0; //restricted static
  Vec3_t x_pred, vr_pred;
  Vec3_t ref_tg;
	#pragma omp parallel for schedule (static) private(P1,P2,end,vr,dist, delta_tg, delta_,delta, x_pred, imp_force, ref_tg, vr_pred, du, normal_cf, m, inside,i,j,crit,dt_fext,kij,omega,psi_cont,e,tgforce,tgvr,norm_tgvr,tgdir,atg) num_threads(Nproc)
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
		for (int i = 0; i<part_id) {
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
        
        x_pred = Particles[P1]->x + Particles[P1]->v * deltat + Particles[P1]->a * deltat * deltat/2.0;
        vr_pred = Particles[P1]->v + Particles[P1]->a * deltat - Particles[P2]->v;
        
        dist =  dot (Particles[P2]->normal, x_pred ) - trimesh[m]-> plane -> pplane;
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
            kij = 2.0 * Particles[P1]->Mass / (deltat * deltat);
						//kij = PFAC * Particles[P1]-> cont_stiff;
						omega = sqrt (kij/Particles[P1]->Mass);
						psi_cont = 2. * Particles[P1]->Mass * omega * DFAC; // Fraser Eqn 3-158
            
            //normal_cf = 2.0 * Particles[P1]->Mass /(deltat*deltat )*delta;

						omp_set_lock(&Particles[P1]->my_lock);
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
						omp_set_lock(&Particles[P1]->my_lock);
						//Particles[P1] -> a += Particles[P1] -> contforce / Particles[P1] -> Mass; 
            Particles[P1] -> a += Particles[P1] -> contforce / Particles[P1] -> Mass; 
						omp_unset_lock(&Particles[P1]->my_lock);
						//cout << "contforce "<<Particles[P1] -> contforce<<endl;


            
					}// if inside
        } //If distance is less than h
    }//Contact Pairs
	}//Nproc
  //cout << "END CONTACT----------------------"<<endl;
	max_contact_force = sqrt (max_contact_force);
	//min_contact_force = sqrt (min_contact_force);
	//cout << "Inside pairs count: "<<inside_geom<<", Inside time: "<<inside_time<<", statically restricted " << stra_restr<<endl;
	int cont_force_count = 0;

}

  
  
};