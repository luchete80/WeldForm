namespace SPH {
//////////////////////////////////////
// HERE PARTICLE DISTRIBUTION IS RADIAL (DIFFERENT FROM PREVIOUS )
// AND HERE distance betwen particles is not even (inner particles are)
// close to each other
// TODO: ANOTHER DOMAIN WITH LESS PARTICLES AT INNER RADIUS POSITIONs}
// ALPHA IN RADIANS
// THIS IS EVEN IN TG 
void Domain::AddCylSliceLengthEven(int tag, double alpha, double Rxy, double Lz, 
																				double r, double Density, double h) {
	//Util::Stopwatch stopwatch;
	std::cout << "\n--------------Generating particles by CylinderBoxLength with defined length of particles-----------" << std::endl;

	size_t PrePS = Particles.Size();
	double xp,yp;
	size_t i,j;
	double qin = 0.03;
	srand(100);
	
	double Lx, Ly;
	
  std::pair <int,Vec3_t> opp_sym; //Position & ID of particles

	//yp=pos;
	int numypart,numxpart;
	int xinc,yinc;
	
	int id_part=0;
	int ghost_rows = 3;
	
	double z0;
	//if (symlength) 	z0 = r;
	//else						
                    z0 = -Lz/2. + r; //CHECK: -Lz/2. - r or -Lz/2.?
	
  int radcount = (Rxy - sqrt(2.)*r) / (2. * r ); //center particle is at (r,r,r)
  cout << "Radial Particle count " <<radcount<<endl;
  
	int part_per_row=0;
  std::vector <int> symm_x;
  std::vector <int> symm_y;
  int x_ghost_per_row = 0;
  //Cal
  int tgcount = alpha* Rxy /(2. * r);
  double dalpha = alpha / tgcount;
  cout << "Tg Particle count " <<tgcount<<endl;
  
  if (Dimension==3) {
    	//Cubic packing
		double zp;
		size_t k=0;
		zp = z0;
		//Calculate row count for non ghost particles
		while (zp <= (z0 + Lz -r)){
			k++; zp += 2.0 * r;			
		}
		//cout << "Particle Row count: "<< k << endl;
		int last_nonghostrow = k;
		k = 0;zp = z0;
    cout << "Length particle count "<<last_nonghostrow+1<<endl;
    
		while (zp <= ( z0 + Lz - r)) {
      for (int alphai = 0; alphai < tgcount + 1; alphai++ ){
        for (int ri = 0; ri < radcount + 1;ri++){
            xp = sqrt(2.)*r + ri * r * cos (alphai*dalpha);
            yp = sqrt(2.)*r + ri * r * sin (alphai*dalpha);
            int tgcount = alpha * (ri+1)*r / (2.*r); 
              
            Particles.Push(new Particle(tag,Vec3_t(xp,yp,zp),Vec3_t(0,0,0),0.0,Density,h,false));
            

         }
         
      } //alpha
			k++;
			zp += 2.0 * r;
		}
    
    cout << "Particles per row: "<<part_per_row<<endl;
    cout << " symmetric particles x, y :"<< symm_x.size()<<", "<<symm_y.size()<<endl;
		//TODO: CONVERT THIS IN A QUARTER CYL SECTOR FUNCTION 
		//Is it convenient to allocate these particles at the end? 		
		//Allocate Symmetry particles, begining from x, y and z
		cout << "Creating ghost particles"<<endl;
    
    
    ///// X AND Y PLANE //////////
    ///// HERE INSERT EACH 2 DIFFERENT SET OF X AND Y GHOST PARTICLES
    // AT THE SAME TIME (GHOST Y PLANE (X DIR INCREMENT) STRAIGHT AND GHOST X PLANE (Y INCREMENT) IS TRANSPOSED)
		zp = z0; k= 0;
		//cout << "zmax"<<( z0 + Lz - r)<<endl;
    int sym_y_count = 0;
    int sym_x_count;

		double Vol = M_PI * Rxy * Rxy * Lz;		
		//double Mass = Vol * Density / (Particles.Size()-PrePS);
		double Mass = Vol * Density /Particles.Size();
		
		cout << "Total Particle count: " << Particles.Size() <<endl;
		cout << "Particle mass: " << Mass <<endl;

		#pragma omp parallel for num_threads(Nproc)
		#ifdef __GNUC__
		for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
		#else
		for (int i=0; i<Particles.Size(); i++)//Like in Domain::Move
		#endif
		{
			Particles[i]->Mass = Mass;
		}

	}//Dim 3

	R = r;									
}


void Domain::AddCylSliceLength(int tag, double alpha, double Rxy, double Lz, 
																				double r, double Density, double h) {
	//Util::Stopwatch stopwatch;
	std::cout << "\n--------------Generating particles by CylinderBoxLength with defined length of particles-----------" << std::endl;

	size_t PrePS = Particles.Size();
	double xp,yp;
	size_t i,j;
	double qin = 0.03;
	srand(100);
	
	double Lx, Ly;
	
  std::pair <int,Vec3_t> opp_sym; //Position & ID of particles

	//yp=pos;
	int numypart,numxpart;
	int xinc,yinc;
	
	int id_part=0;
	int ghost_rows = 2;
	
	double z0;
	//if (symlength) 	z0 = r;
	//else						
                    z0 = -Lz/2. + r; //CHECK: -Lz/2. - r or -Lz/2.?
	
  int radcount = (Rxy - sqrt(2.)*r) / (2. * r ); //center particle is at (r,r,r)
  cout << "Radial Particle count " <<radcount<<endl;
  
	int part_per_row=0;
  std::vector <int> symm_x;
  std::vector <int> symm_y;
  int x_ghost_per_row = 0;
  //Cal
  int tgcount;

  cout << "Tg Particle count " <<tgcount<<endl;
  
  int part_count = 0;
  
  if (Dimension==3) {
    	//Cubic packing
		double zp;
		size_t k=0;
		zp = z0;
		//Calculate row count for non ghost particles
		while (zp <= (z0 + Lz -r)){
			k++; zp += 2.0 * r;			
		}
		//cout << "Particle Row count: "<< k << endl;
		int last_nonghostrow = k;
		k = 0;zp = z0;
    cout << "Length particle count "<<last_nonghostrow+1<<endl;

  //tangential particles
  int plane_ghost_part_1[last_nonghostrow][radcount][ghost_rows]; //First plane id
  int plane_ghost_part_2[last_nonghostrow][radcount][ghost_rows]; //First plane id
  
  //plane_ghost_part_1_count[last_nonghostrow][radcount]; // Not always is possible to count for ghost count
  
  //ref particles
  //int plane_ghost_part_3[][]
  
    //First increment is in radius
		while (zp <= ( z0 + Lz - r)) {
      int rcount = 0; //Used only for ghost count
      for (double ri = 0. ; ri < Rxy; ri += 2.*r){
        double rtot = ri + sqrt(2.) * r;
        //cout << "ri "<<ri<<endl;
        
        double dalpha;
        if (ri == 0.) {tgcount =1; dalpha = 0.;}
        else {

          tgcount = (ceil)((alpha* ri )/(2. * r)) + 1;  
          dalpha = alpha / (tgcount-1);         
          //cout << "tg count "<<tgcount<<", dalpha"<<dalpha<<", alpha ri"<<alpha * ri<<"ri "<<ri <<endl;
        }
        for (int alphai = 0; alphai < tgcount; alphai++ ){
            xp = sqrt(2.) * r + ri * cos (alphai*dalpha);
            yp = sqrt(2.) * r + ri * sin (alphai*dalpha);
            Particles.Push(new Particle(tag,Vec3_t(xp,yp,zp),Vec3_t(0,0,0),0.0,Density,h,false));
            if (alphai < ghost_rows){
              plane_ghost_part_1[k][rcount][alphai] = part_count++; //plane_ghost_part_1[last_nonghostrow][radcount][ghost_rows]
              plane_ghost_part_2[k][rcount][alphai] = part_count++; //plane_ghost_part_1[last_nonghostrow][radcount][ghost_rows]
            }
            part_count++;
            
            // if ( i < ghost_rows ){ //X PLANE SYMMETRY
              // symm_x.push_back(id_part);
              // if (k==0) x_ghost_per_row++;
              // //if (k==0) 
              // //  Particles[id_part]->ID = id_part; //ONLY FOR TESTING IN PARAVIEW!
            // }
            // if ( j < ghost_rows) { //Y PLANE SYMMETRY
              // symm_y.push_back(id_part);
              // //if (k==0) 
              // //  Particles[id_part]->ID = id_part; //ONLY FOR TESTING IN PARAVIEW!
            // }
            // if (zp == z0)
              // part_per_row++;
            
            // id_part++;
            // xp += 2.*r;
          // }
          // yp += 2.*r;
          // yinc +=1;

         }
        rcount++;
      } //alpha
			k++;
			zp += 2.0 * r;
		}
    
    cout << "Particles per row: "<<part_per_row<<endl;
    cout << " symmetric particles x, y :"<< symm_x.size()<<", "<<symm_y.size()<<endl;
		//TODO: CONVERT THIS IN A QUARTER CYL SECTOR FUNCTION 
		//Is it convenient to allocate these particles at the end? 		
		//Allocate Symmetry particles, begining from x, y and z
		cout << "Creating ghost particles"<<endl;
    
    
    ///// X AND Y PLANE //////////
    ///// HERE INSERT EACH 2 DIFFERENT SET OF X AND Y GHOST PARTICLES
    // AT THE SAME TIME (GHOST Y PLANE (X DIR INCREMENT) STRAIGHT AND GHOST X PLANE (Y INCREMENT) IS TRANSPOSED)
		zp = z0; k= 0;
		//cout << "zmax"<<( z0 + Lz - r)<<endl;
    int sym_y_count = 0;
    int sym_x_count;
		// while (zp <= ( z0 + Lz - r)) {

			// numypart = numpartxy;	//And then diminish by 2 on each y increment
			// yinc = 1;	//particle row from the axis
			// //cout << "y particles: "<<numypart<<endl;
			// for (j=0; j < ghost_rows ; j++){//DY INCREMENT 
				// xp = r;
				// yp = - r - 2*r*(yinc -1); //First increment is radius, following ones are 2r			
				// numxpart = calcHalfPartCount(r, Rxy, yinc);
				// //cout << "x particles: "<<numypart<<endl;
        // sym_x_count = x_ghost_per_row *k + j;
        // //THE INCREMENTS ARE IN X, i.e. PARALLEL TO Y PLANE
				// for (i=0; i < numxpart;i++) { //DX INCREMENT
					// //if (random) Particles.Push(new Particle(tag,Vec3_t((x + qin*r*double(rand())/RAND_MAX),(y+ qin*r*double(rand())/RAND_MAX),(z+ qin*r*double(rand())/RAND_MAX)),Vec3_t(0,0,0),0.0,Density,h,Fixed));
					// //	else    
					// Particles.Push(new Particle(tag,Vec3_t(xp,yp,zp),Vec3_t(0,0,0),0.0,Density,h,false)); //First insert on y plane 
					// Particles.Push(new Particle(tag,Vec3_t(yp,xp,zp),Vec3_t(0,0,0),0.0,Density,h,false)); //Transpose for x plane

          // //if (k==0) 
          // //  Particles[id_part  ]->ID = symm_y[sym_y_count]; //ONLY FOR TESTING IN PARAVIEW!  
          // //SYMMETRY ON Y PLANE IS ALTERNATED ON INDICES
          // //if (k==0) 
          // //  Particles[id_part+1]->ID = symm_x[sym_x_count]; //ONLY FOR TESTING IN PARAVIEW!           
          
					// Particles[id_part  ]->inner_mirr_part = symm_y[sym_y_count];
					// Particles[id_part+1]->inner_mirr_part = symm_x[sym_x_count];
          
          // GhostPairs.Push(std::make_pair(symm_y[sym_y_count],id_part  ));
          // GhostPairs.Push(std::make_pair(symm_x[sym_x_count],id_part+1));
          
          // Particles[id_part  ]->ghost_plane_axis = 1;
          // Particles[id_part+1]->ghost_plane_axis = 0;
          // Particles[id_part  ]->is_ghost = true;
          // Particles[id_part+1]->is_ghost = true;
          
          // Particles[id_part  ]->not_write_surf_ID = true; //TO NOT BE WRITTEN BY OUTER SURFACE CALC
          // Particles[id_part+1]->not_write_surf_ID = true; //TO NOT BE WRITTEN BY OUTER SURFACE CALC
          // //ONLY FOR TESTING SYMMETRY PLANES!
          // //Particles[id_part  ]->ID = 1; //ONLY FOR TESTING IN PARAVIEW!  
          // //Particles[id_part+1]->ID = 0; //ONLY FOR TESTING IN PARAVIEW!   
					
          // id_part+=2;
					// xp += 2.*r;
          // sym_y_count++;
          // sym_x_count+=ghost_rows;
				// }//x rows
				// yinc++;
			// }//y rows
			// k++;
			// zp += 2.0 * r;	
			// //cout << "zp " <<zp<<endl;
		// }
		
		double Vol = M_PI * Rxy * Rxy * Lz;		
		//double Mass = Vol * Density / (Particles.Size()-PrePS);
		double Mass = Vol * Density /Particles.Size();
		
		cout << "Total Particle count: " << Particles.Size() <<endl;
		cout << "Particle mass: " << Mass <<endl;

		#pragma omp parallel for num_threads(Nproc)
		#ifdef __GNUC__
		for (size_t i=0; i<Particles.Size(); i++)	//Like in Domain::Move
		#else
		for (int i=0; i<Particles.Size(); i++)//Like in Domain::Move
		#endif
		{
			Particles[i]->Mass = Mass;
		}

	}//Dim 3

	R = r;									
}

};