#include "matvec.h" 
#define VMIN_FOR_FRICTION	1.e-3

#define VMAX_FOR_STA_FRICTION	1.e-3
#include "Plane.h"

namespace SPH {

//////////////////////////////// 
//// From Zhan: A SPH framework for dynamic interaction between soil and
////            rigid body system with hybrid contact method
////////////////////////////////
inline void Domain::CalcContactForcesZhan(){
	
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
	Vec3_t stick_cf;
  double normal_cf;
	double delta_ = 0.;
	double deltat_cont;
  double crit;
  
  double kij, omega,psi_cont;
  int i,j;  //For inside testing
	
	int P1,P2;
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
  Vec3_t vp;
  
  int max_reached_part = 0; //TEST
  int sta_frict_particles = 0;
  int stra_restr = 0; //restricted static
	#pragma omp parallel for schedule (static) private(P1,P2,end,vr,vp,delta_,delta, deltat_cont, m, inside,i,j,crit,stick_cf,normal_cf,dt_fext,kij,omega,psi_cont,e,tgvr,norm_tgvr,tgdir,atg) num_threads(Nproc)
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
						// delta = (deltat - deltat_cont) * delta_;
						// DEBUG THINGS, REMOVE ////////////
						inside_part[P1] ++;
            vp = Particles[P1]->v + Particles[P1] -> a * deltat; //Projected body vel if no contact
            //ATTENTION: ASSUMING CONSTANT VELOCITY OF RIGID BODY
            //This results of equalling projected velocity differences to zero
            // NORMAL FORCE
            normal_cf = dot (vp - Particles[P2]->v,	Particles[P2]->normal) * Particles[P1]->Mass / (deltat);
						omp_set_lock(&Particles[P1]->my_lock);
						Particles[P1] -> contforce = normal_cf *  Particles[P2]->normal;  // //Eqn 30. Zhan
						omp_unset_lock(&Particles[P1]->my_lock);

            stick_cf = (vp - Particles[P2]->v) * Particles[P1]->Mass / (deltat) - Particles[P1] -> contforce; //Eqn 31 Zhan
            
            if (norm (stick_cf) < friction * normal_cf){
              omp_set_lock(&Particles[P1]->my_lock);
                Particles[P1] -> a += stick_cf;  // NORMAL DIRECTION, Fraser 3-159
              omp_unset_lock(&Particles[P1]->my_lock);              
            } else { //Sliding
              tgdir = stick_cf/norm (stick_cf);
              omp_set_lock(&Particles[P1]->my_lock);
                Particles[P1] -> a += friction * Particles[P1] -> contforce;  // NORMAL DIRECTION, Fraser 3-159
              omp_unset_lock(&Particles[P1]->my_lock);   
            }

						// if   (force2 > max_contact_force ) max_contact_force = force2;
						// else if (force2 < min_contact_force ) min_contact_force = force2;
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


}; //SPH
