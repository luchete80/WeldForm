#include "matvec.h" 
#define VMIN_FOR_FRICTION	1.e-3

#define VMAX_FOR_STA_FRICTION	1.e-3
#include "Plane.h"


namespace SPH {

//////////////////////////////// 
//// From Wang: Simulating frictional contact in smoothed particle hydrodynamics
//// https://link.springer.com/article/10.1007/s11431-013-5262-x
//// Wang, Wu, GU, HUA, Science China 2013
////////////////////////////////
inline void Domain::CalcContactForcesLS(){

	double min_force_ts_=1000.;
// https://stackoverflow.com/questions/10850155/whats-the-difference-between-static-and-dynamic-schedule-in-openmp

	////////////////////////
	// DEBUG THINGS //////
	std::vector <int> inside_part(Particles.Size());
	double min_delta,max_delta;
	min_delta = 1000.; max_delta = 0.;
	int inside_time,inside_geom;

  #pragma omp parallel for num_threads(Nproc)
  for (int i = 0;i<Particles.Size();i++){
		//omp_set_lock(&Particles[i]->my_lock);
		Particles[i] -> contforce = 0.; //RESET    
		Particles[i] -> delta_cont = 0.; //RESET    
		Particles[i] -> tgdir = 0.;				//TODO: DELETE (DEBUG) 
    Particles[i] -> q_fric_work = 0.;
    Particles[i] -> cshearabs = 0.; //cshear module
		//omp_unset_lock(&Particles[i]->my_lock);
		inside_part[i] = 0;
		inside_time=inside_geom=0;
    Particles[i] -> q_cont_conv = 0.;
    Particles[i] -> friction_hfl = 0.;
  }

  for (int m=0;m<meshcount;m++) tot_cont_heat_cond[m] = 0.;
  
  
	max_contact_force = 0.;
	double min_contact_force = 1000.;
	int inside_pairs = 0;
	double delta_ = 0.;
  double delta;
  double crit;
  
  double kij, omega,psi_cont;
  int i,j;  //For inside testing
	
	int P1,P2;
  Vec3_t tgforce, tgforce_dyn;
  Vec3_t imp_force;

  Vec3_t vr;
  double dt_fext;
  Element* e;
  bool inside;
  
  Vec3_t tgvr, tgdir;
  double norm_tgvr;
  Vec3_t delta_tg;
  double max_vr = 0.;
  int m;
  double abs_fv;
  
  Vec3_t du;// If no contact
 
  ext_forces_work_step = 0.;
  //USE TRUE: FALSE DOES NOT WORK
  //bool ref_accel = false;   //true: tg force is compared to current tg accel

  double fr_sta, fr_dyn;
            
  Vec3_t atg;
  bool end;
  contact_force_sum = 0.;
  contact_reaction_sum = 0.;
  double dist;
  int max_reached_part = 0; //TEST
  int sta_frict_particles = 0;
  int stra_restr = 0; //restricted static
  Vec3_t x_pred, vr_pred;
  Vec3_t ref_tg;
  double dS2;

	#pragma omp parallel for schedule (static) private(P1,P2,end,vr,dist, delta_tg, delta_,delta, abs_fv, x_pred, imp_force, fr_sta, fr_dyn, ref_tg, vr_pred, du, m, inside,i,j,crit,dt_fext,kij,omega,psi_cont,e,tgforce,tgforce_dyn,tgvr,norm_tgvr,tgdir,atg, dS2) num_threads(Nproc)
  //tgforce
	#ifdef __GNUC__
	for (size_t k=0; k<Nproc;k++) 
	#else
	for (int k=0; k<Nproc;k++) 
	#endif	
	{
		Vec3_t xij, qj;
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
        
        x_pred = Particles[P1]->x + Particles[P1]->v * deltat + Particles[P1]->a * deltat * deltat/2.0;
        vr_pred = Particles[P1]->v + Particles[P1]->a * deltat - Particles[P2]->v;
        
        dist =  dot (Particles[P2]->normal, x_pred ) - trimesh[m]-> element[Particles[P2]->element] -> pplane;
        //cout << "dist "<<dist<<", h "<<Particles[P1]->h<<endl;
        if( dist  < Particles[P1]->h) {

          qj = Particles[P1]->x - dist * Particles[P2]->normal;
                                 //Check if it is inside triangular element
					//Find a vector 
					//Fraser 3-147
					inside = true;
          if (trimesh[m]->dimension == 3){
            i=0;		
            while (i<3 && inside){
              j = i+1;	if (j>2) j = 0;
              crit = dot (cross ( *trimesh[m]->node[e -> node[j]] 
                                          - *trimesh[m]->node[e -> node[i]],
                                          qj  - *trimesh[m]->node[e -> node[i]]),
                                Particles[P2]->normal);
              if (crit < 0.0) inside = false;
              i++;
            }
          } else { //MESH DIMENSION = 2
            i=0;
            while (i<2 && inside){
              j = i+1;	if (j>1) j = 0;
              crit = dot ( *trimesh[m]->node[e -> node[j]] 
                                          - *trimesh[m]->node[e -> node[i]],
                                          qj  - *trimesh[m]->node[e -> node[i]]);
              if (crit < 0.0) inside = false;
              i++;
            }
          }
					
					if (inside ) { //Contact point inside element, contact proceeds
            inside_geom++;
            end=true;

            delta = Particles[P1]->h - dist;
            //cout << "dist "<<dist<<", h "<<Particles[P1]->h<< ", delta "<<delta<<endl;
						// DAMPING
						//Calculate SPH and FEM elements stiffness (series)
						//Since FEM is assumed as rigid, stiffness is simply the SPH one 
            // if (!gradKernelCorr)
              kij = 2.0 * Particles[P1]->Mass / (deltat * deltat) * PFAC;
            // else //TESTING PHASE
              // kij = Particles[P1]->Mass / (deltat * deltat);
						
            //kij = PFAC * Particles[P1]-> cont_stiff;
						omega = sqrt (kij/Particles[P1]->Mass);
						psi_cont = 2. * Particles[P1]->Mass * omega * DFAC; // Fraser Eqn 3-158

						omp_set_lock(&Particles[P1]->my_lock);
              Particles[P1] -> contforce = (kij * delta - psi_cont * delta_) * Particles[P2]->normal; // NORMAL DIRECTION, Fraser 3-159    
              //Particles[P1] -> delta_cont = delta;
						omp_unset_lock(&Particles[P1]->my_lock);
            
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
                omp_set_lock(&dom_lock);
								this -> min_force_ts = min_force_ts_;
                omp_unset_lock(&dom_lock);
						}
						// omp_set_lock(&Particles[P1]->my_lock);
            // Particles[P1] -> a += Particles[P1] -> contforce / Particles[P1] -> Mass; 
						// omp_unset_lock(&Particles[P1]->my_lock);

            
            fr_sta = friction_sta;
            fr_dyn = friction_dyn;
            if (friction_function == Linear)
              fr_sta = fr_dyn = friction_m * Particles[P1] ->T + friction_b;
            
            Vec3_t fT(0.,0.,0.);
            
            if (fr_sta > 0.) { 
              Vec3_t fricold = Particles[P1] -> tgforce;
              
              Vec3_t kdeltae = PFAC * Particles[P1]->Mass*vr/deltat;
              
              //double fy = friction_mu*glm::length(fN);	//coulomb friction
              double fy = fr_sta * norm(Particles[P1] -> contforce);  
              Vec3_t fstar = fricold - kdeltae;

              if (norm(fstar) > fy) {
                fT  = fy*fstar/norm(fstar);
              } else {
                fT = fstar;
              }
              
              //Particles[P1] -> tgforce = fT;
              
            omp_set_lock(&Particles[P1]->my_lock);
              Particles[P1] -> contforce += fT; // NORMAL DIRECTION, Fraser 3-159    
						omp_unset_lock(&Particles[P1]->my_lock);

						// omp_set_lock(&Particles[P1]->my_lock);
						// //Particles[P1] -> a += Particles[P1] -> contforce / Particles[P1] -> Mass; 
            // Particles[P1] -> a += Particles[P1] -> tgforce / Particles[P1] -> Mass; 
						// omp_unset_lock(&Particles[P1]->my_lock);  


            }//friction
            
            omp_set_lock(&dom_lock);            
              contact_force_sum += norm(Particles[P1] ->contforce);
              //contact_force_sum_v += Particles[P1] ->contforce;
              contact_reaction_sum += dot (Particles[P1] -> a,Particles[P2]->normal)* Particles[P1]->Mass;
              //ext_forces_work_step += dot (Particles[P1] -> contforce,//Particles[P2]->v);
              //ext_forces_work_step += /*dot(*/Particles[P1] -> contforce/*,1./norm(Particles[P2]->v)*Particles[P2]->v)*/ * Particles[P2]->v; //Assuming v2 and forces are parallel
              ext_forces_work_step += dot(Particles[P1] -> contforce,Particles[P2]->v);
            omp_unset_lock(&dom_lock);	

            
					}// if inside
        } //If distance is less than h
    }//Contact Pairs
    
    //contact_force_sum = norm (contact_force_sum_v);
	}//Nproc
  //cout << "END CONTACT----------------------"<<endl;
	max_contact_force = sqrt (max_contact_force);
	//min_contact_force = sqrt (min_contact_force);
	//cout << "Inside pairs count: "<<inside_geom<<", Inside time: "<<inside_time<<", statically restricted " << stra_restr<<endl;
	int cont_force_count = 0;
	
  ext_forces_work += ext_forces_work_step * deltat;
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
  for (int m=0;m<meshcount;m++) 
    accum_cont_heat_cond +=tot_cont_heat_cond[m] * deltat;
}


}; //SPH
