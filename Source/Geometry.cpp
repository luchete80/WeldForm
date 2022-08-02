#include "Plane.h"
#include <vector>

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

    //reflected particles
    //Particles.Push(new Particle(tag,-Particles[0]->x,Vec3_t(0,0,0),0.0,Density,h,false));
    Vec3_t xr;
    for (int k=0;k<last_nonghostrow; k++){
      for (int ri=0;ri < ghost_rows; ri++){
        for (int alphai = 0; alphai < tgcount_ref_ghost[ri]; alphai++ ){
          int id_part = plane_ghost_part_3[k][ri][alphai]; 
          xr = - Particles[id_part]->x;
          xr(2) = Particles[id_part]->x(2);
          Particles.Push(new Particle(tag,xr,Vec3_t(0,0,0),0.0,Density,h,false));
          GhostPairs.Push(std::make_pair(id_part,part_count));
          //ONLY FOR DEBUG
          // if (k==last_nonghostrow-1)
            // Particles[part_count  ]->ID = id_part;    
          //Plane Ghost is not necessary here
          Particles[part_count  ]-> ghost_type = Mirror_XY;
          Particles[part_count  ]-> is_ghost = true;
          Particles[part_count  ]->not_write_surf_ID = true;   
              
          part_count++;
          ghost_count++;
        }
      }
    }
          
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

};