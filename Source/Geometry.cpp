#include "Plane.h"
#include <vector>

#include "NastranVolReader.h"

namespace SPH {
//////////////////////////////////////
// HERE PARTICLE DISTRIBUTION IS RADIAL (DIFFERENT FROM PREVIOUS )
// AND HERE distance betwen particles is not even (inner particles are)
// close to each other

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
    z0 = /*-Lz/2. + */ r; //CHECK: -Lz/2. - r or -Lz/2.?
	
  int radcount = Rxy / (2. * r ); //center particle is at (r,r,r)
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
  int plane_ghost_part_1[2][last_nonghostrow][radcount][ghost_rows]; //First plane id
  //plane_ghost_part_count[3][radcount][]; // Not always is possible to count for ghost count
  
  //Reflected particles (in relation to center), same angle on -x,-y coordinates
  int plane_ghost_part_3[last_nonghostrow][ghost_rows + 1 ] [ghost_rows + 1 ];
  
  for (int i=0;i<2;i++)
  for (int k=0;k<last_nonghostrow;k++)
  for (int r=0;r<radcount;r++)
  for (int a=0;a<ghost_rows;a++)
    plane_ghost_part_1[i][k][r][a] = -1;
   
  int tgcount_ref_ghost[ghost_rows];
  
    //First increment is in radius
		while (zp <= ( z0 + Lz - r)) {
      int rcount = 0; //Used only for ghost count
      for (double ri = 0. ; ri < Rxy; ri += 2.*r){
        //cout << "ri "<<ri<<endl;
        
        double dalpha;
        if (ri == 0.) {tgcount =1; dalpha = 0.;}
        else {

          tgcount = (ceil)((alpha* ri )/(2. * r)) + 1;  
          dalpha = alpha / (tgcount-1);         
          //cout << "tg count "<<tgcount<<", dalpha"<<dalpha<<", alpha ri"<<alpha * ri<<"ri "<<ri <<endl;
          if (rcount > 0 && rcount < ghost_rows +1) tgcount_ref_ghost[rcount - 1] = tgcount;
        }
        for (int alphai = 0; alphai < tgcount; alphai++ ){
            xp =  /*r +*/ ri * cos (alphai*dalpha);
            yp =  /*r +*/ ri * sin (alphai*dalpha);
            Particles.Push(new Particle(tag,Vec3_t(xp,yp,zp),Vec3_t(0,0,0),0.0,Density,h,false));
            if (ri == 0.){
              Particles[part_count]->is_fixed = true;
            }
            //A particle can be on all zones 
            if (alphai < ghost_rows){
              // if (k==5)
              // cout << "pl0 k, r, alpha, id_part " << k<<", " <<rcount <<", " <<alphai<<", " <<part_count<<endl;
              plane_ghost_part_1[0][k][rcount][alphai] = part_count; //plane_ghost_part_1[last_nonghostrow][radcount][ghost_rows]
            }
            if ((tgcount - alphai - 1) < ghost_rows){ // Opposite plane
              // if (k==5)
                // cout << "pl1 k, r, alpha, id_part " << k<< ", " <<rcount <<", " <<tgcount - alphai - 1<<", " <<part_count<<endl;
              plane_ghost_part_1[1][k][rcount][tgcount - alphai - 1] = part_count; //plane_ghost_part_1[last_nonghostrow][radcount][ghost_rows]
            } //Reflected particles
            if (rcount > 0 && rcount < ghost_rows +1){ //Particle at origin x,y are not reflected
              plane_ghost_part_3[k][rcount-1][alphai] = part_count;
            }
            Particles[part_count  ]->not_write_surf_ID = true;   
            //Only for debug  
            // if (k==last_nonghostrow-1){
              // Particles[part_count  ] ->ID = part_count; 
              // Particles[part_count  ]->not_write_surf_ID = true;   
            // }//ONLY FOR DEBUG
            
            part_count++;
         }
        rcount++;
      } //alpha
			k++;
			zp += 2.0 * r;
		}
    cout << "Not Ghost Particle Count: "<<part_count<<endl;
    cout << "Particles per row: "<<part_per_row<<endl;
		cout << "Creating ghost particles"<<endl;
    
    int ghost_count = 0;
    
    ///// X AND Y PLANE //////////
    ///// HERE INSERT EACH 2 DIFFERENT SET OF X AND Y GHOST PARTICLES
    // AT THE SAME TIME (GHOST Y PLANE (X DIR INCREMENT) STRAIGHT AND GHOST X PLANE (Y INCREMENT) IS TRANSPOSED)
		zp = z0; k= 0;
		//cout << "zmax"<<( z0 + Lz - r)<<endl;
    //OUTER, NOT COORDINATE NORMAL
    Vec3_t normal_2[2], tg[2][2]; //tg is [plane][comp]
    double pplane[2];
    
    //Planes are not exactly at 0,0, they have an offset of r of normal dir
    normal_2 [0] = Vec3_t(0., -1.0 ,0.);
    pplane [0]   = dot (normal_2[0],Particles[0]->x + normal_2[0] * r);
    tg[0][0] = Vec3_t(1., 0. ,0.);
    tg[0][1] = tg[1][1] = Vec3_t(0.,0.,1.);    
    tg[1][0] = Vec3_t(cos(alpha), sin(alpha),0.);
    
    normal_2 [1] = Vec3_t(cos(alpha + M_PI/2.), sin(alpha + M_PI/2.),0.);
    pplane [1]   = dot (normal_2[1],Particles[0]->x + normal_2[1] * r);
    
    for (int i=0;i<2;i++){
      //Plane *p = new Plane(normal_2[i],tg[i][0],tg[i][1],pplane[i]);
      //Plane *p = new Plane();
      planes.push_back(new Plane(normal_2[i],tg[i][0],tg[i][1],pplane[i]));
    }
    
    cout << "pplane"<<pplane<<endl;

    for (int k=0;k<last_nonghostrow;k++){
      for (int ri=0;ri<radcount;ri++){
        for (int alphai = 0; alphai < ghost_rows; alphai++ ){
          for (int i=0;i<2;i++){
            int id_part = plane_ghost_part_1[i][k][ri][alphai]; 
            if (alphai == 0){ //Only First row 
              Particles[id_part  ]   ->plane_ghost = planes[i];
              Particles[id_part  ]   ->correct_vel_acc = true;
            }
          }
        }
      }
    }
    
    int sym_y_count = 0;
    int sym_x_count;
    for (int k=0;k<last_nonghostrow;k++){
      for (int ri=0;ri<radcount;ri++){
        for (int alphai = 0; alphai < ghost_rows; alphai++ ){
          for (int i=0;i<2;i++){
            int id_part = plane_ghost_part_1[i][k][ri][alphai]; 
            if (id_part != -1){
              //cout << "k, r, alpha, id part << " << k<<"," <<ri<< ", " <<alphai << ", " <<id_part<<endl;
              Particles[id_part  ]->ghost_plane_axis = 1;
              
              //Move particle across normal
              //Distance to plane: t = (n . X) - pplane
              //Being pplane the plane coeff n . X = pplane              
              double dist = dot (Particles[id_part]->x,normal_2[i] ) - pplane[i] ;
              Vec3_t xg = Particles[id_part]->x + 2 * normal_2[i] * abs(dist); //Outer normal
              xp = xg(0); yp = xg(1); 
              
              zp = Particles[id_part]->x(2);
            //if (i==1 && alphai == 0) { //TEST
              Particles.Push(new Particle(tag,Vec3_t(xp,yp,zp),Vec3_t(0,0,0),0.0,Density,h,false));
              GhostPairs.Push(std::make_pair(id_part,part_count));
              Particles[part_count  ]->is_ghost = true;
              Particles[part_count  ]->plane_ghost = planes[i];
              Particles[part_count  ]->ghost_type = Symmetric;
              //Only for debug
              // if (k==last_nonghostrow-1)
                // Particles[part_count  ]->ID = id_part;         
                    
              Particles[part_count  ]->not_write_surf_ID = true;               
              part_count++;
              ghost_count++;
            } //If !=-1
            //}//TEST
          }//i          
        }//alpha i
      }//r
    }//z
    
    //
    //reflected particles
    // Vec3_t xr;
    // for (int k=0;k<last_nonghostrow; k++){
      // for (int ri=0;ri < ghost_rows; ri++){
        // for (int alphai = 0; alphai < tgcount_ref_ghost[ri]; alphai++ ){
          // int id_part = plane_ghost_part_3[k][ri][alphai]; 
          // xr = - Particles[id_part]->x;
          // xr(2) = Particles[id_part]->x(2);
          // Particles.Push(new Particle(tag,xr,Vec3_t(0,0,0),0.0,Density,h,false));
          // GhostPairs.Push(std::make_pair(id_part,part_count));
          // //ONLY FOR DEBUG
          // // if (k==last_nonghostrow-1)
            // // Particles[part_count  ]->ID = id_part;    
          // //Plane Ghost is not necessary here
          // Particles[part_count  ]-> ghost_type = Mirror_XY;
          // Particles[part_count  ]-> is_ghost = true;
          // Particles[part_count  ]->not_write_surf_ID = true;   
              
          // part_count++;
          // ghost_count++;
        // }
      // }
    // }
          
    cout << "Ghost count "<<ghost_count<<endl;
    
		
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



void Domain::AddCylUniformLength(int tag, double Rxy, double Lz, 
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
    z0 = /*-Lz/2. + */ r; //CHECK: -Lz/2. - r or -Lz/2.?
	
  int radcount = Rxy / (2. * r ); //center particle is at (r,r,r)
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

    //First increment is in radius
		while (zp <= ( z0 + Lz - r)) {
      int rcount = 0; //Used only for ghost count
      for (double ri = 0. ; ri < Rxy; ri += 2.*r){
        //cout << "ri "<<ri<<endl;
        
        double dalpha;
        if (ri == 0.) {tgcount =1; dalpha = 0.;}
        else {

          tgcount = (ceil)((2.*M_PI* ri )/(2. * r));  
          dalpha = 2.*M_PI / (tgcount);         

        }
        for (int alphai = 0; alphai < tgcount; alphai++ ){
            xp =  /*r +*/ ri * cos (alphai*dalpha);
            yp =  /*r +*/ ri * sin (alphai*dalpha);
            Particles.Push(new Particle(tag,Vec3_t(xp,yp,zp),Vec3_t(0,0,0),0.0,Density,h,false));
            
            part_count++;
         }
        rcount++;
      } //alpha
			k++;
			zp += 2.0 * r;
		}

		double Vol = M_PI * Rxy * Rxy * Lz;		
		//double Mass = Vol * Density / (Particles.Size()-PrePS);
		double Mass = Vol * Density /Particles.Size();
		
		cout << "Total Particle count: " << Particles.Size() <<endl;
		cout << "Particle mass: " << Mass <<endl;
    
    solid_part_count = Particles.Size();

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

//In case elements were generated at Nastran Node Positions
// void Domain::GenerateSPHMesh(const int &tag, NastranVolReader &nr,double Density){

  // for (int n=0;n<nr.node_count;n++){
    // Vec3_t v(nr.node[3*n],nr.node[3*n+1],nr.node[3*n+2]);
    // double h=1.;
    // Particles.Push(new Particle(tag,v,Vec3_t(0,0,0),0.0,Density,h,false));
  // }
   // cout << "Generated "<<Particles.Size()<<" particles. "<<endl;
   
  // std::vector<double> h(Particles.Size());
  // std::vector<double> elxnod(Particles.Size());
  
  // // //cout << "Normals"<<endl;
  // for (int e=0;e<nr.elem_count;e++){ //Look for element
    // Vec3_t v = nr.elem[e]->centroid;

    // double dist[8];
    // //
    // for (int i=0;i < nr.elem[e]->node_count; i++){
      // //Vec3_t v[3];
      // int k = i+1;
      // if (i==2) k = 0; 
      // dist[i] = norm(Particles[nr.elcon[3*e + k]]->x - Particles[nr.elcon[3*e+i]]->x );
      // elxnod[nr.elcon[3*e + k]]++;
      // elxnod[nr.elcon[3*e + i]]++;
      // h[nr.elcon[3*e + k]] += dist[i];
      // h[nr.elcon[3*e + i]] += dist[i];
    // }
  // }  
  // for (int i=0;i<Particles.Size();i++)
    // Particles[i]->h = h[i]/elxnod[i];

    // Particles.Push(new Particle(tag,v,Vec3_t(0,0,0),0.0,Density,h,false));
  
// }

//Generated at centroids
//ORIGINAL, MORE RAM
void Domain::GenerateSPHMesh(const int &tag, NastranVolReader &nr,double Density, double hfac){

   for (int e=0;e < nr.elem.size();e++){
    //Vec3_t v(nr.node[3*n],nr.node[3*n+1],nr.node[3*n+2]);
    double h=1.;
    //cout << "centroid "<<nr.elem[e]->centroid<<endl; 
    Particles.Push(new Particle(tag,nr.elem[e]->centroid,Vec3_t(0,0,0),0.0,Density,h,false));
  }
   cout << "Generated "<<Particles.Size()<<" particles. "<<endl;
   
  std::vector<double> h(Particles.Size());
  std::vector<int> el_nod_count(nr.node_count);
  short el_nod[nr.node_count][50]; //element list per node, TODO: CHANGE BY ELEMENT POINTER
  std::vector<int> el_nb_count(Particles.Size());
  //el_nod.resize();
  
  // //cout << "Normals"<<endl;
  //First check each element node and add current element as shared 
  int maxn = 0;
  for (int e=0;e<nr.elem_count;e++){ //Look for element
    for (int i=0;i < nr.elem[e]->nodecount; i++){
      int nid = nr.elem[e]->node[i];
      el_nod[nid][el_nod_count[nid]]= e;
      
      // if (nid == 30) cout << nr.elem[e]->id << "ncount " <<i<< " internal el id " <<e<<endl;
      // if (nid == 30) cout << "el_nod[nid][el_nod_count[nid]]" <<el_nod[nid][el_nod_count[nid]]<<endl;
      el_nod_count[nid] ++;
      if (e==142 || nr.elem[e]->id == 7679)cout << nid << " ";
      if (el_nod_count[nid] > maxn)
        maxn = el_nod_count[nid];
    }
  }
  cout << "Max elements shared by nodes: "<<maxn<<endl;
  cout << "Nopde 30 "<<endl;
  
  // for (int i=0;i<el_nod_count[30];i++)
    // cout << nr.elem[el_nod[30][i]]->id<< " " ;
  // cout <<endl;
  
  // for (int n = 0 ; n<100;n++){
    // cout << el_nod_count[n] << ": "<<endl;
    // for (int nn=0;nn<el_nod_count[n];nn++)
      // cout << nr.elem[el_nod[n][nn]]->id << " ";
    // cout << endl;
  // }
  
  // for (int n=0;n<10;n++)
    // cout << "
  //Now check from the list which elements are shared by node
  maxn = 0;
  int emax;
  for (int e=0;e<nr.elem_count;e++){ //Look for element
    for (int i=0;i < nr.elem[e]->nodecount; i++){//Element node
      int nid = nr.elem[e]->node[i];
      for (int en = 0; en < el_nod_count[nid]; en++){//Each element shared this node
        if(el_nod[nid][en] != e){
          //el_nb[e][el_nb_count[e]] = el_nod[nid][en];
          el_nb_count[e]++;
          if (el_nb_count[e]++ > maxn){
            maxn = el_nb_count[e];
            emax = e;
          }
          h[e] += norm (nr.elem[el_nod[nid][en]]->centroid - nr.elem[e]->centroid); 
          if (e==142){
            cout << "node " <<i << ", id "<<nid <<", el " <<nr.elem[el_nod[nid][en]]->id<<endl;
            
          }
        }
      }//Loop through shared element nodes 
    }
  }
  //Loop throug nodes 
  cout << "Max elem neighbours " <<maxn << " for elem id" << nr.elem[emax]->id <<endl;

  
  double hmax = 0.;
  int imax, imin;
  double hmin =1000.;
  for (int i=0;i<Particles.Size();i++){
    Particles[i]->h = h[i]/el_nb_count[i] * hfac;
    Particles[i]->Nb = el_nb_count[i];
    if (Particles[i]->h > hmax){
      hmax = Particles[i]->h;
      imax = i;
    }
    else if (Particles[i]->h < hmin){
      hmin = Particles[i]->h;
      imin = i;
    }
  }
  cout << "h max is " << hmax << ", in particle " << imax << ", elem id " << nr.elem[imax]->id<<endl;
  cout << "h min is " << hmin << ", in particle " << imin << ", elem id " << nr.elem[imin]->id<<endl;
  cout << "nb count " <<  el_nb_count[imax]<<endl;
  // for (int nb=0;nb<el_nb_count[imax];nb++){
    // cout << nr.elem[el_nb[imax][nb]]->id<<", ";
  // }
  cout <<endl;
  
  double volmax = 0.;
  double tot_mass = 0.;
  double tot_vol = 0.;
  
  for (int i=0; i < Particles.Size(); i++){
    Particles[i]->Mass = Density * nr.elem[i]->vol;
    tot_mass += Particles[i]->Mass;
    tot_vol += nr.elem[i]->vol;;
    if (nr.elem[i]->vol > volmax )
      volmax = nr.elem[i]->vol;
  }
  cout << "Max Element volume "<< volmax << endl;
  cout << "Total Volume: " << tot_vol;
  cout << "Total Mass: " << tot_mass << endl;
  
  solid_part_count = Particles.Size();
  
}

// //Slow,but does not crash, it uses less RAM
// void Domain::GenerateSPHMesh(const int &tag, NastranVolReader &nr,double Density, double hfac){

  // for (int e=0;e < nr.elem.size();e++){
    // //Vec3_t v(nr.node[3*n],nr.node[3*n+1],nr.node[3*n+2]);
    // double h=1.;
    // //cout << "centroid "<<nr.elem[e]->centroid<<endl;
    // Particles.Push(new Particle(tag,nr.elem[e]->centroid,Vec3_t(0,0,0),0.0,Density,h,false));
  // }
   // cout << "Generated "<<Particles.Size()<<" particles. "<<endl;
   
  // std::vector<double> h(Particles.Size());
  // std::vector<int> el_nod_count(nr.node_count);
  // short el_nod[nr.node_count][50]; //element list per node, TODO: CHANGE BY ELEMENT POINTER
  // std::vector<int> el_nb_count(Particles.Size());
  // //el_nod.resize();
  
  // // //cout << "Normals"<<endl;
  // //First check each element node and add current element as shared 
  // int maxn = 0;
  // int emax;
  // for (int e=0;e<nr.elem_count;e++){ //Look for element
    // for (int i=0;i < nr.elem[e]->nodecount; i++){
      // int nid = nr.elem[e]->node[i];
      // for (int en = 0 ; en < nr.elem_count;en++) { //Loop though all elements looking for this node
        // for (int in=0;in < nr.elem[en]->nodecount; in++){
          // if (nr.elem[e]->node[i] == nr.elem[en]->node[in] ) {
            // el_nb_count[e]++;
          // }
          // if (el_nb_count[e]++ > maxn){
            // maxn = el_nb_count[e];
            // emax = e;
          // }
        // }
      // }
    // }
  // }
   // cout << "Max elem neighbours " <<maxn << " for elem id" << nr.elem[emax]->id <<endl;
  
  // double hmax = 0.;
  // int imax;
  // for (int i=0;i<Particles.Size();i++){
    // Particles[i]->h = h[i]/el_nb_count[i];
    // Particles[i]->Nb = el_nb_count[i];
    // if (Particles[i]->h > hmax){
      // hmax = Particles[i]->h;
      // imax = i;
    // }
  // }
  // cout << "h max is " << hmax << ", in particle " << imax << ", elem id " << nr.elem[imax]->id<<endl;
  // cout << "nb count " <<  el_nb_count[imax]<<endl;
  // // for (int nb=0;nb<el_nb_count[imax];nb++){
    // // cout << nr.elem[el_nb[imax][nb]]->id<<", ";
  // // }
  // cout <<endl;
  
  // for (int i=0;i<Particles.Size();i++){
    // Particles[i]->Mass = Density * nr.elem[i]->vol;
  // }
  
// }


};