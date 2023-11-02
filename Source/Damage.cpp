using namespace std;

//S. Chakraborty, A. Shaw / International Journal of Impact Engineering 58 (2013) 84e95

namespace SPH {

///// calculate dam_D factor for each pair
///// 
inline void Domain::CalcDamage(){
  // Particle *PP[2];
  // //int pp[2];
    
  // double eps_max;
  // Mat3_t sig;
  
	// #pragma omp parallel for schedule (static) private (PP) num_threads(Nproc)
	// #ifdef __GNUC__
	// for (size_t k=0; k<Nproc;k++) 
	// #else
	// for (int k=0; k<Nproc;k++) 
	// #endif	
	// {
  // for (size_t p=0; p<SMPairs[k].Size();p++) {
    
    // // //#ifndef NONLOCK_SUM
    // // if (!nonlock_sum){
    // PP[0]	= Particles[SMPairs[k][p].first];
    // PP[1]	= Particles[SMPairs[k][p].second];	

    // // pt[0]	= SMPairs[k][p].first;
    // // pt[1]	= SMPairs[k][p].second;	
    
    // Vec3_t xij	= PP[0]->x - PP[1]->x;
		// Vec3_t vij	= PP[0]->v - PP[1]->v;
    
    // double rij	= norm(xij);
    
    // // CHECK INITIATION CRITERIA
    // if (dam_D[k][p] == 0.0){ //OR dam_df0[][] == 0.0
      // //Average ef_max for each particle? 
      // if ((rij - dam_r0_unif)/dam_r0_unif >= eps_max){
        // dam_rf0[k][p] = rij;
      // }//INIT CRITERION
    // } else {
			// if (dam_D[k][p] == 0.0 < 1) {
			// //EVOLUTION
				// double siRR[] = {0.0,0.0}; //Fom each particle
				// double delta = rij - dam_rf0[k][p];
				// Vec3_t xij2 = xij * xij;
				// double f = 1.0/(rij*rij);
				// double urij = -f * (dot(xij,vij));
				// double ci[2];
				// double ffi[2];
				// for (int i=0;i<2 ;i++){
					// //P1->Sigma eqn. 24
					// sig = PP[i]->Sigma;
					// for (int j=0;j<3;j++) siRR[i] += xij2(j)*sig(j,j);
					
					// siRR[i] += 2.0*(sig(0,1)*xij(0)*xij(1) + sig(1,2)*xij(1)*xij(2) + sig(2,0)*xij(2)*xij(0));
					// siRR[i] *= f; 
					// ci[i] = SoundSpeed(PP[i]->PresEq, PP[i]->Cs, PP[i]->Density, PP[i]->RefDensity); //TODO: CALCULATE THIS ONCE! (IN THE ACCEL CALC IS ALREADY DONE)
					// ffi[i] = PP[i]->Density * ci[i];
				// }//i
				// double sijrr = ( siRR[0] * ffi[1] + siRR[1] * ffi[0] ) / (ffi[0] + ffi[1]); //Eqn 22.
				
				// //In the paper, dam_rf0 is delta = 0|t_0 in the Fig. 3
				// //Points of the cohesive curve are (dam_rf0, sigma_max) and (deltamax, 0)
				// //ASSUME p1 and p2 are the same material because they are united....
				// double m = -PP[0]->mat->damage->sigma_max / (PP[0]->mat->damage->delta_max - dam_rf0[k][p]); //SLOPE
				// //mx+b = 0 (m*deltamax  + b = 0)
				// double b = -m * PP[0]->mat->damage->delta_max;
				// //Calculate max transferable traction from linear cohesive law
				// double sig_al_t = m * delta + b;
				// if (sijrr > sig_al_t) {
					// if (delta < PP[0]->mat->damage->delta_max)
						// dam_D[k][p] = 1.0 - sig_al_t / sijrr;
					// else dam_D[k][p] = 1.0;
				// }
				
			// } //if D < 1
    // }//Evolution
    // // }//nonlock sum -> only way coded
  // }//for pair p
  // } //proc k
}

}; //SPH