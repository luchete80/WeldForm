/***********************************************************************************
* PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids        * 
*             and soils) using Smoothed Particle Hydrodynamics method              *   
* Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo-Torres              *
*                                                                                  *
* This file is part of PersianSPH                                                  *
*                                                                                  *
* This is free software; you can redistribute it and/or modify it under the        *
* terms of the GNU General Public License as published by the Free Software        *
* Foundation; either version 3 of the License, or (at your option) any later       *
* version.                                                                         *
*                                                                                  *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY  *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.         *
*                                                                                  *
* You should have received a copy of the GNU General Public License along with     *
* PersianSPH; if not, see <http://www.gnu.org/licenses/>                           *
************************************************************************************/

#include "Domain.h"
#include "Input.h"

#include "Input.h"
#include "InteractionAlt.cpp"
#include "Mesh.h"

#include "SolverFraser.cpp"
#include "Geometry.cpp"
#include "SolverKickDrift.cpp"
#include "SolverLeapfrog.cpp"

#define TAU		0.005
#define VMAX	10.0

#define PRINTVEC(v)	cout << v[0]<<", "<<v[1]<<", "<<v[2]<<endl;

using namespace std;
using namespace SPH;

// template <typename T > 
// void ReadParameter(T in, nlohmann::json in) {
  // cout << "Reading Configuration parameters..."<<endl; 
  // cout << "Time step size: ";
  // readValue(in["timeStepSize"], /*scene.timeStepSize*/ts);
  // cout << ts << endl;
// }



void UserAcc(SPH::Domain & domi)
{
	double vcompress;

	if (domi.getTime() < TAU ) 
		vcompress = VMAX/TAU * domi.getTime();
	else
		vcompress = VMAX;
	
	#pragma omp parallel for schedule (static) num_threads(domi.Nproc)
	#ifdef __GNUC__
	for (size_t i=0; i<domi.Particles.Size(); i++)
	#else
	for (int i=0; i<domi.Particles.Size(); i++)
	#endif
	{
		for (int bc=0;bc<domi.bConds.size();bc++){
			if (domi.Particles[i]->ID == domi.bConds[bc].zoneId ) {
				if (domi.bConds[bc].type == Velocity_BC ){ //VELOCITY
					if (domi.bConds[bc].valueType == 0) {
            domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
            domi.Particles[i]->v		= domi.bConds[bc].value;          
          } else if (domi.bConds[bc].valueType == 1) {///amplitude
            for (int j=0;j<domi.amps.size();j++){
              if(domi.amps[j].id == domi.bConds[bc].ampId){
                double val = domi.bConds[bc].ampFactor * domi.amps[j].getValAtTime(domi.getTime());
                Vec3_t vec = val * domi.bConds[bc].value;
                domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
                domi.Particles[i]->v		= vec;
              }//if if match
            }//for amps
          } //VALUE TYPE == AMPLITUDE 
				}//  == VELOCITY
        else if (domi.bConds[bc].type == Temperature_BC){
          domi.Particles[i]->T = domi.bConds[bc].T;
        } else if (domi.bConds[bc].type == Convection_BC){
          
        }
        
			}//ZoneID 			
		}//BC
	}
  
  if (domi.contact){
    for (int bc=0;bc<domi.bConds.size();bc++){
      for (int m=0;m<domi.trimesh.size();m++){
        if (domi.trimesh[m]->id == domi.bConds[bc].zoneId){
          if (domi.bConds[bc].valueType == 0) { ///constant
            domi.trimesh[m]->SetVel(domi.bConds[bc].value);
            domi.trimesh[m]->SetRotAxisVel(domi.bConds[bc].value_ang);
          }
          else if (domi.bConds[bc].valueType == 1) {///amplitude
            for (int i=0;i<domi.amps.size();i++)
              if(domi.amps[i].id == domi.bConds[bc].ampId){
                double val = domi.bConds[bc].ampFactor * domi.amps[i].getValAtTime(domi.getTime());
                Vec3_t vec = val * domi.bConds[bc].value;
                // cout << "Time, vec, ampfactor"<<domi.getTime()<< ", "<<vec<<", amp "<<domi.bConds[bc].ampFactor<<endl;
                // cout << "bc val "<<domi.bConds[bc].value<<endl;
                // cout << "amp val " <<domi.amps[i].getValAtTime(domi.getTime())<<endl;
                domi.trimesh[m]->SetVel(vec);
              }
 				// readValue(bc["amplitudeId"], 		bcon.ampId);
				// readValue(bc["amplitudeFactor"], 	bcon.ampFactor);           
          }
        }//zoneID
      }//mesh
    }//bcs
    
  }//contact
}

size_t findLastOccurrence(string str, char ch)
{
 
    // To store the index of the result
    size_t found;

    found = str.rfind(ch);
    // If string doesn't have
    // character ch present in it
    if (found == string::npos) {
        cout << "Character " << ch
             << " is not present in"
             << " the given string.";
    }
 
    // Else print the position
    else {
        cout << "The last occurrence of '"
             << ch << "' is found at index: "
             << found << endl;
    }
  return found;
}
 

int main(int argc, char **argv) try {

	if (argc > 1) {
		string inputFileName=argv[1];	
		std::ifstream i(argv[1]);
		json j;
		i >> j;

		nlohmann::json config 		= j["Configuration"];
		nlohmann::json material 	= j["Materials"];
		nlohmann::json domblock 	= j["DomainBlocks"];
		nlohmann::json domzones 	= j["DomainZones"];
		nlohmann::json amplitudes 	= j["Amplitudes"];
		nlohmann::json rigbodies 		= j["RigidBodies"];
    nlohmann::json contact_ 		= j["Contact"];
		nlohmann::json bcs 			= j["BoundaryConditions"];
		nlohmann::json ics 			= j["InitialConditions"];

		
		SPH::Domain	dom;
		dom.filename = inputFileName;
    size_t pos = 0;
    size_t test = findLastOccurrence(inputFileName, '\\');
    if (test != string::npos) pos = test;
    //cout << "pos of json "<<inputFileName.find(".json")<<endl;
    string out_name = inputFileName.substr(pos, inputFileName.find(".json") - pos +1) + "out";
    //cout << "Out file: "<< out_name << endl;
    dom.out_file.open(out_name.c_str()/*, std::ofstream::out | std::ofstream::app*/);
		dom.Dimension	= 3;
		
		string kernel;
    double ts;
    
    ostringstream oss_out;
    
    oss_out << "-----------------------------------------------"<<endl;
    oss_out << "------------------ WELDFORM -------------------"<<endl;
    oss_out << "----------------- v. 0.4.2.4 ------------------"<<endl;
    oss_out << "------------------ 20241127 -------------------"<<endl;
    oss_out << "-- 331c3a8ed9250581180e967fcb8845fa8c8aaa91 ---"<<endl;
    oss_out << "-----------------------------------------------"<<endl<<endl<<endl;
    
    cout << oss_out.str();
    dom.out_file << oss_out.str();
    
    cout << "Reading Configuration parameters..."<<endl; 
    
    int np = 4;
    readValue(config["Nproc"], np);   
    dom.Nproc	= np;
        
    string sumType = "Nishimura";
    readValue(config["sumType"], /*scene.timeStepSize*/sumType);
    if (sumType == "Locking") dom.nonlock_sum = false;
    else if (sumType == "Nishimura") dom.nonlock_sum = true;
    else cout << "sumType value not valid. Options are \"Locking\" and \"Nishimura\". "<<endl; 
    
    
		cout << "Time step size: ";
    readValue(config["timeStepSize"], /*scene.timeStepSize*/ts);
    cout << ts << endl;
    dom.Kernel_Set(Qubic_Spline);
    
    
    string solver = "Mech";
    readValue(config["solver"],solver);
    
    if (solver=="Thermal" ||solver=="Mech-Thermal" || solver=="Mech-Thermal-KickDrift" || solver=="Mech-Thermal-Fraser" || solver=="Mech-Thermal-LeapFrog")
      dom.thermal_solver = true;
		
		readValue(config["integrationMethod"], dom.Scheme); //0:Verlet, 1:LeapFrog, 2: Modified Verlet

     	//dom.XSPH	= 0.5; //Very important

        double dx,r,h;


		readValue(config["particleRadius"], r);
		double hfactor;
		readValue(config["hFactor"], hfactor);
    bool h_update = false;
    dom.GeneralAfter = & UserAcc;
		
	
		//////////////
		// MATERIAL //
		//////////////
		double rho,E,nu,K,G,Cs,Fy;
    Fy = -1.0;
    double Et, Ep;  //Hardening (only for bilinear and multilear)
    std::vector<double> c;
    c.resize(10);
    std::vector<double> e_range (2,0.0);
    std::vector<double> er_range(2,0.0);
    std::vector<double> T_range (2,0.0);
    e_range[1]=er_range[1]=T_range[1]=1.0e10;
    
    string mattype = "Bilinear";
    bool plastic_heat = false;
    cout << "Reading Material.."<<endl;
    cout << "Type.."<< endl; readValue(material[0]["type"], 		mattype);
    cout << "Density.."<< endl; readValue(material[0]["density0"], 		rho);
    readValue(material[0]["youngsModulus"], 	E);
    readValue(material[0]["poissonsRatio"], 	nu);
    readValue(material[0]["yieldStress0"], 	Fy);
    readArray(material[0]["const"], 		c);
    readArray(material[0]["strRange"],  e_range );
    readArray(material[0]["strdotRange"], er_range);
    readArray(material[0]["tempRange"],   T_range );
    readValue(material[0]["plasticHeat"], plastic_heat);
    
    dom.pl_heating = plastic_heat;
    
    Material_ *mat; //SINCE DAMAGE MATERIAL, MATERIAL ALWAYS HAS TO BE ASSIGNED
    Elastic_ el(E,nu);
    cout << "Mat type  "<<mattype<<endl;
    if      (mattype == "Bilinear")    {
      Ep = E*c[0]/(E-c[0]);		                              //only constant is tangent modulus
      cout << "Material Constants, Et: "<<c[0]<<endl;
			mat = new Material_(el);	//THIS IS NEW WITH DAMAGE; 
														//AND IS COHERENT WITH BILINEAR MATERIAL WHICH DID NOT HAVE ITS OWN Material Class    
		} else if (mattype == "Hollomon")    {
      mat = new Hollomon(el,Fy,c[0],c[1]);
      cout << "Material Constants, K: "<<c[0]<<", n: "<<c[1]<<endl;
    } else if (mattype == "JohnsonCook") {
      //Order of input is: [A,B,n,C,eps_0,m,Tm,Tt] //FIRST STRAIN, THEN STRAIN RATE AND THEN THERMAL
      /////INPUT IN CONSTRUCTOR IS A,B,C,
                               //A ,B,,n, c,eps_0,m,T_m, T_transition
      mat = new JohnsonCook(el,c[0],c[1],c[2],c[3],c[4], c[5],c[6],c[7]); //First is hardening // A,B,C,m,n_,eps_0,T_m, T_t);	 //FIRST IS n_ than m
      cout << "Material Constants, A: "<<c[0]<<", B: "<<c[1]<<", n: "<<c[2]<<"C: "<<c[3]<<", eps_0: "<<c[4]<<"m: "<<c[5]<<", T_m: "<<c[6]<<", T_t: "<<c[7]<<endl;
    } else if (mattype == "GMT") {
      //Order of input is: n1,n2  c1,c2, m1,m2, I1, I2
      mat = new GMT(el,c[0],c[1],c[2],c[3],c[4], c[5],c[6],c[7],
                       e_range [0],e_range [1],
                       er_range[0],er_range[1],
                       T_range [0],T_range [1]); //First is hardening // A,B,C,m,n_,eps_0,T_m, T_t);	 //FIRST IS n_ than m
      cout << "GMT Material Constants: "<<endl<<
                                  "n1: "<<c[0]<<", n2: "<<c[1]<<endl<<
                                  "c1: "<<c[2]<<", c2: "<<c[3]<<endl<<
                                  "m1: "<<c[4]<<", m2: "<<c[5]<<endl<<
                                  "I1: "<<c[6]<<", I2: "<<c[7]<<endl;
    }    
    
    else                              throw new Fatal("Invalid material type.");
    
    
		string damage_mod="";
		readValue(material[0]["damageModel"],damage_mod);
		DamageModel *damage;
    if      (damage_mod == "Rankine"){
			double smax, Gf;
			readValue(material[0]["smax"], 				smax);
			readValue(material[0]["fracEnergy"], 	Gf);
			damage= new RankineDamage(smax,Gf);
			mat->damage = damage;
			cout << "Assigned Rankine Damage Model"<<endl;
			dom.model_damage = true;
			dom.nonlock_sum = false;
		} else if      (damage_mod == "JohnsonCook"){
			
			double smax, Gf;
      std::vector <double> D(5);
			readArray(material[0]["damageParams"], 		D);
			damage= new JohnsonCookDamage(D[0],D[1],D[2],D[3],D[4],mat->getRefStrainRate()); //Correct this
			mat->damage = damage;
			cout << "Assigned Johnson Cook Damage Model"<<endl;
			cout << "Damage Model Params: "<<D[0]<<", "<<D[1]<<", "<<D[2]<<", "<<D[3]<<", "<<D[4]<<endl;
			dom.model_damage = true;
			dom.nonlock_sum = false;
		}
    if (damage_mod==""){
      cout << "NO DAMAGE MODEL SET."<<endl;
    } else if(damage_mod=="Rankine" || damage_mod=="JohnsonCook"){
      damage->mat = mat;
    }
			
    // THERMAL PROPERTIES

    double k_T, cp_T, th_ex;
    th_ex = 0.0;
    readValue(material[0]["thermalCond"], 	  k_T);
    readValue(material[0]["thermalHeatCap"], 	cp_T);    
    readValue(material[0]["thermalExp"], 	  th_ex); //Expansion
    
    cout << "Thermal Parameters: "<<endl;
    cout << "Expansion: "<<th_ex<<endl;
    cout << "HeatCap:" <<cp_T<<endl;
    cout << "thermalCond"<<k_T<<endl;
    
    cout << "Done. "<<endl;
       
		K= E / ( 3.*(1.-2*nu) );
		G= E / (2.* (1.+nu));

		dx 	= 2.*r;
    h	= dx*hfactor; //Very important
    Cs	= sqrt(K/rho);

        double timestep,cflFactor;
		int cflMethod;
    double output_time;
    double sim_time;
    string cont_alg = "Wang";
    bool auto_ts[] = {true, false, false}; //IN SOME METAL CUTS HAS BEEN FOUND THAT ONLY VELOCITY CRITERIA DIVERGES.
		readValue(config["cflMethod"], cflMethod);

		if (cflMethod == 0)
			readValue(config["timeStepSize"], timestep);
		else {
			readValue(config["cflFactor"], cflFactor);
			timestep = (cflFactor*h/(Cs));
		}
    double min_ts_manual = timestep/1.e-5;
    readValue(config["minTS"], min_ts_manual);
      dom.manual_min_ts = min_ts_manual;
    readValue(config["outTime"], output_time);
    readValue(config["simTime"], sim_time);
    readBoolVector(config["autoTS"], auto_ts);
    double alpha = 1.;
    double beta = 0.;
    bool h_upd = false;
    double tensins = 0.3;
    int nb_upd_freq = 5;
    bool kernel_grad_corr = false;
    int gradType = 0;
    readValue(config["artifViscAlpha"],alpha);
    readValue(config["artifViscBeta"],beta);
    readValue(config["contAlgorithm"],cont_alg);
    readValue(config["kernelGradCorr"],kernel_grad_corr);
    readValue(config["stressGradType"],gradType);
    readValue(config["smoothlenUpdate"],h_upd);
    readValue(config["nbsearchFreq"],nb_upd_freq);
    dom.auto_ts = auto_ts[0];
    dom.auto_ts_acc = auto_ts[1];
    dom.auto_ts_cont = auto_ts[2];
    
    Gradient_Type gt = gradType;
    dom.Gradient_Approach_Set(gradType);
    //dom.GradientType = gradType;
		
    readValue(config["tensileInstability"],tensins);
    if (h_upd) //default is false...
      dom.h_update = true;
    /////////////-/////////////////////////////////////////////////////////////////////////////////
		// DOMAIN //
		////////////
		Vec3_t start,L;
    int id;
		string domtype = "Box";
    int matID;
    string gridCS = "Cartesian";
    double slice_ang = 2.0000001 * M_PI;
    bool sym[] = {false,false,false};
		readValue(domblock[0]["id"], 	id);
		readVector(domblock[0]["start"], 	start);
		cout << "Reading Domain dim" << endl;  readVector(domblock[0]["dim"], 	L);
		cout << "Reading Domain type" << endl; readValue(domblock[0]["type"], 	domtype); //0: Box
    cout << "Reading Domain mat id" << endl;  readValue(domblock[0]["matID"], 	matID); //0: Box
    cout << "Grid Coordinate System" << endl;  readValue(domblock[0]["gridCoordSys"], 	gridCS); //0: Box
    cout << "Slice Angle " << endl;  readValue(domblock[0]["sliceAngle"], 	slice_ang); //0: Box
    readBoolVector(domblock[0]["sym"], 	sym); //0: Box
    for (int i=0;i<3;i++){ //TODO: Increment by Start Vector
			dom.DomMax(0) = L[i];
			dom.DomMin(0) = -L[i];
		}		


		
		// inline void Domain::AddCylinderLength(int tag, Vec3_t const & V, double Rxy, double Lz, 
									// double r, double Density, double h, bool Fixed) {
												
		//dom.AddCylinderLength(1, Vec3_t(0.,0.,-L/10.), R, L + 2.*L/10. + dx, r, rho, h, false); 
		
    if (abs(L[2]) < h ) {
      dom.Dimension = 2;
      cout << "Z Value is less than h. Dimension is set to 2. "<<endl;
      cout << "Dimension also could be set in config section." <<endl;
    }
    
		cout << "Dimensions: "<<endl;
		PRINTVEC(L)
		if (domtype == "Box"){
      cout << "Adding Box ..."<<endl;  
      if ( gridCS == "AxiSymmetric") {
        dom.dom_bid_type = AxiSymmetric;
        cout << "PROBLEM TYPE: 2D AXISYMMETRIC"<<endl;
      }
			dom.AddBoxLength(id ,start, L[0] , L[1],  L[2] , r ,rho, h, 1 , 0 , false, false );		
			cout << "Solid Part count "<<dom.solid_part_count<<endl;

		}
		else if (domtype == "Cylinder"){
      cout << "Adding Cylinder";      
			if (sym[0] && sym[1]){
        cout << " with symmetry..."<<endl;
        dom.AddXYSymCylinderLength(0, L[0]/2., L[2], r, rho, h, false, sym[2]); 
      }
      else {
        cout << "..."<<endl;
        if ( gridCS == "Cartesian"){
          cout << "PROBLEM TYPE: 3D"<<endl;
          cout << "DIM: "<<dom.Dimension<<endl;
          dom.AddCylinderLength(0, start, L[0]/2., L[2], r, rho, h, false, sym[2]); 
        }
        else if (gridCS == "Cylindrical"){
          if (slice_ang==2.0 * M_PI){
            dom.AddCylUniformLength(0, L[0]/2.,L[2], r, rho, h);
          } else {
            dom.AddCylUniformLength(0, L[0]/2.,L[2], r, rho, h, M_PI/16.0 , 1, L[0]/4.0); 
            dom.dom_bid_type = AxiSymm_3D;

          }
        } else if (gridCS == "CylRadial"){
          if (slice_ang==2.0 * M_PI){
            //dom.AddCylSliceLength(0, L[0]/2.,L[2], r, rho, h);
            cout << "ERROR. RADIAL CYL SHOULD BE LESS THAN 2PI"<<endl;
          } else {
            dom.AddCylSliceLength(0, L[0]/2.,L[2], r, rho, h, slice_ang /*, 1, L[0]/4.0*/); 
            dom.dom_bid_type = AxiSymm_3D;

          }          
        }
      }
    } else if (domtype == "File") {
        double scalefactor = 1.0d;
        readValue(domblock[0]["scaleFactor"],scalefactor);
        string filename = "";
        readValue(domblock[0]["fileName"], 	filename); 
        cout << "Reading Particles Input file " << filename <<endl;  
        dom.ReadFromLSdyna(filename.c_str(), rho);
        
        cout << "Scaling by factor: "<< scalefactor<<endl;
        for (int i=0;i<dom.Particles.Size();i++)
            dom.Particles[i]->x *= scalefactor;
        
        Vec3_t translation;
        readVector(domblock[0]["translation"], 	      translation);      //Or value linear
        cout << "Translation vector: "<<translation<<endl;
        for (int i=0;i<dom.Particles.Size();i++)
            dom.Particles[i]->x += translation;        
            
        double totmass = 0.0;
        readValue(domblock[0]["totMass"], 	totmass); 
        
        bool calcParticleRadius = false;
        readValue(domblock[0]["calcParticleRadius"], 	calcParticleRadius); 
        
        Vec3_t dims = dom.getBboxDims();
        
        cout << "--------DOMAIN DIMENSION SET TO: ";
        if (dims[2]>0.0){
          dom.Dimension = 3;
        } else {
          dom.Dimension = 2;          
        }
        cout << dom.Dimension<<endl;
        
        
        if (calcParticleRadius){
          cout << "calculating avg distance ..."<<endl;
          double avgdist = dom.getAvgMinDist();
          cout << "Avg particle distance "<<avgdist<<endl;
        
          cout <<"Setting smoothing length to "<<avgdist<<endl;
          for (int i=0;i<dom.Particles.Size();i++){
              dom.Particles[i]->h = avgdist*hfactor;
          }
        } else {
          cout << "Using default Smoothing Length of: "<<h<<endl;
          for (int i=0;i<dom.Particles.Size();i++){
              dom.Particles[i]->h = h;
          }          
          
        }
        if (totmass != 0){
        double mass = totmass/dom.Particles.Size();
        cout << "Appliyng particle mass: "<<mass<<endl;
        for (int i=0;i<dom.Particles.Size();i++)
            dom.Particles[i]->Mass = mass;

        } else {
          
          cout << "ERROR. TOT  MASS UNKNOWN"<<endl;
        }
        
        
          //x =dom_d->x_h[p].x;
          //y =dom_d->x_h[p].y;
          //z = dom_d->x_h[p].z;
          //dom.Particles.push_back(new SPH::Particle(0,Vector(x,y,z),Vector(0,0,0),0.0,rho,h,false));
          //dom.Particles[p]->Mass = dom_d->m_h[p];
         
          //if (dom_d->realloc_ID)dom.Particles[p]->ID = dom_d->ID_h[p];
          //tot_mass+=dom_d->m_h[p];
        //}
        //delete dom_d->x_h,dom_d->m_h;
        //printf( "Total Mass Readed from LS-Dyna: %fn", tot_mass);      
      
      
    }//File

        cout <<"t  			= "<<timestep<<endl;
        cout <<"Cs 			= "<<Cs<<endl;
        cout <<"K  			= "<<E<<endl;
        cout <<"G  			= "<<nu<<endl;
        cout <<"Fy 			= "<<Fy<<endl;
		cout <<"dx 			= "<<dx<<endl;
		cout <<"h  			= "<<h<<endl;
		cout <<"-------------------------"<<endl;
		cout <<	"Dim: "<<dom.Dimension<<endl;				
		cout << "Particle count: "<<dom.Particles.Size()<<endl;



    cout << "Domain Zones "<<domzones.size()<<endl;		
		for (auto& zone : domzones) { //TODO: CHECK IF DIFFERENTS ZONES OVERLAP
			// MaterialData* data = new MaterialData();
			int zoneid;
			Vec3_t start,end;
			readValue(zone["id"], 		zoneid);
			readVector(zone["start"], 	start);
			readVector(zone["end"], 	end);
      cout << "Zone id "<<zoneid<<endl;
			// cout << "Dimensions: "<<endl;
			// PRINTVEC(start)
			// PRINTVEC(end)
			int partcount =dom.AssignZone(start,end,zoneid);
      std::cout<< "Zone "<<zoneid<< ", particle count: "<<partcount<<std::	endl;
		}
    
    //////////////////////////////////////////////////////////
    ////////////////// RIGID BODIES //////////////////////////
    std::vector<TriMesh *> mesh;
    int mesh_count = 0;
    for (int rb = 0; rb < rigbodies.size();rb++){
      string rigbody_type;
      bool contact = false;
      if (readValue(rigbodies[rb]["type"],rigbody_type))
        contact = true;
      
      if (contact){
        Vec3_t dim;
      
                
        readVector(rigbodies[rb]["start"], 	start);       
        readVector(rigbodies[rb]["dim"], 	dim); 
        bool flipnormals = false;
        readValue(rigbodies[rb]["flipNormals"],flipnormals);
        
        double heatcap = 1.;
        readValue(rigbodies[rb]["thermalHeatCap"],heatcap);
        cout << "Reading Contact surface "<<endl;
        //TODO: WRitE TO PArTiclES
        if (rigbody_type == "File"){
          // string filename = "";
          // readValue(rigbodies[0]["fileName"], 	filename); 
          // cout << "Reading Mesh input file..." << endl;
          // SPH::NastranReader reader("Tool.nas", flipnormals);
        }
        else {
          if (dim (0)!=0. && dim(1) != 0. && dim(2) !=0. && rigbody_type == "Plane")
            throw new Fatal("ERROR: Contact Plane Surface should have one null dimension");
        }
        
        cout << "Set contact to ";

        cout << "true."<<endl;
        dom.contact = true;
        cout << "Reading contact mesh..."<<endl;
        //TODO: CHANGE TO EVERY DIRECTION
        int dens = 10;
        readValue(rigbodies[rb]["partSide"],dens);
        if (rigbody_type == "Plane"){
          // TODO: CHECK IF MESH IS NOT DEFINED
          mesh.push_back(new TriMesh); 
          mesh[mesh_count]->AxisPlaneMesh(2, flipnormals, start, Vec3_t(start(0)+dim(0),start(1)+dim(1), start(2)),dens);
        } else if (rigbody_type == "Line"){
          mesh.push_back(new TriMesh);
          mesh[mesh_count]->dimension = 2;
          //if (dim(0)>0.0 && dim(0)>1.0) cout << "ERROR. 
          if (dim(0)>0.0)
            mesh[mesh_count]->AxisPlaneMesh(1, flipnormals, start, Vec3_t(start(0)+dim(0),start(1)+dim(1), 0.0),dens);        
          else if (dim(1)>0.0)
            mesh[mesh_count]->AxisPlaneMesh(0, flipnormals, start, Vec3_t(start(0)+dim(0),start(1)+dim(1), 0.0),dens); 
          else 
            cout << "ERROR. Line has null dimension."<<endl;
          cout << "Rigid Body start pos: "<<start(0)+dim(0)<<", "<<start(1)+dim(1)<<endl;
        } else if (rigbody_type == "File"){
          
          string filename = "";
          readValue(rigbodies[rb]["fileName"], 	filename); 
          cout << "Reading Mesh input file " << filename <<endl;
          SPH::NastranReader reader(filename.c_str());
            mesh.push_back (new SPH::TriMesh(reader,flipnormals ));
        }

        double scalefactor = 1.0d;
        readValue(rigbodies[rb]["scaleFactor"],scalefactor);
        if (scalefactor != 1.0){
          cout << "Scaling mesh..."<<endl;
          mesh[mesh_count]->Scale(scalefactor);
        }

        Vec3_t translation(0,0,0);
        readVector(rigbodies[rb]["translation"], 	      translation);      //Or value linear
        cout << "Moving Mesh by vector "<<translation<<endl;
        if (norm(translation)>1.0e-6){
        //  mesh[mesh_count]->Move(translation);
          for (int n=0;n<mesh[mesh_count]->node.Size();n++){
            *(mesh[mesh_count]->node[n]) += translation;
          }
          mesh[mesh_count]->CalcCentroids();
          mesh[mesh_count]->CalcNormals(); 
          //NOT UPDATING COEFFS SINCE NFAR IS NOT CALC YET (IS CALC BELOW IN SPHERES) 
        }
        //After calc Normals show avg normal pointing:
        Vec3_t avgn = 0.;
        for (int e = 0; e < mesh[mesh_count]->element.Size(); e++) 
          avgn += mesh[mesh_count]->element[e] -> normal;
        avgn /=mesh[mesh_count]->element.Size();
        cout << "Surface avg Normal Direction (check to see if flipNormals needs to be set): "<<endl<<avgn<<endl;
        
        cout << "Creating Spheres.."<<endl;
        //mesh.v = Vec3_t(0.,0.,);
        mesh[mesh_count]->CalcSpheres(); //DONE ONCE
        double hfac = 1.1;	//Used only for Neighbour search radius cutoff
        cout << "Adding mesh particles ...";
        int id;
        readValue(rigbodies[rb]["zoneId"],id);
        cout << "Id: "<<id <<endl;
        dom.AddTrimeshParticles(mesh[mesh_count], hfac, id); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
        dom.m_contact_force.push_back(Vec3_t(0.0,0.0,0.0));
        double penaltyfac = 0.5;
        std::vector<double> fric_sta(1), fric_dyn(1), heat_cond(1);
        readValue(contact_[0]["fricCoeffStatic"], 	fric_sta[0]); 
        readValue(contact_[0]["fricCoeffDynamic"], 	fric_dyn[0]); 
        readValue(contact_[0]["heatCondCoeff"], 	  heat_cond[0]);
        
        readValue(contact_[0]["penaltyFactor"], 	penaltyfac); 
        
        mesh_count ++;

        // readValue(rigbodies[0]["contAlgorithm"],cont_alg);
        cout << "Contact Algorithm set to ";
        if (cont_alg == "Seo") {
          cout << "SEO"<<endl;
          dom.contact_alg = Seo;
        } else if (cont_alg == "LSDyna") {
          dom.contact_alg = LSDyna;
          cout << "LS_Dyna"<<endl;
        } else if (cont_alg == "Fraser") { /////FRASER NOT WORKING
          // dom.contact_alg = Fraser;
          cout << "Wang. ATTENTION: Fraser  algorithm does not working"<<endl;
        }else {
          cout << "Wang"<<endl;
        }
        
        //cout << "Contact Algortihm: "<< cont_alg.c_str() <<end;
        
        bool heat_cond_ = false;
        bool heat_fric_ = false;
        readValue(contact_[0]["heatConductance"], 	heat_cond_);
        readValue(contact_[0]["heatFriction"], 	heat_fric_);
        if (heat_cond_) {
          dom.cont_heat_cond = true;
          dom.contact_hc = heat_cond[0];
        }

        if (heat_fric_) {
          dom.cont_heat_fric = true;
        }
        

        dom.friction_dyn = fric_dyn[0];
        dom.friction_sta = fric_sta[0];
        cout << "Contact Friction Coefficients, Static: "<<dom.friction_sta<<", Dynamic: "<< dom.friction_sta<<endl;
        
        dom.PFAC = penaltyfac;
        dom.DFAC = 0.0;
        cout << "Contact Penalty Factor: "<<dom.PFAC<<", Damping Factor: " << dom.DFAC<<endl; 
        cout << "Contact heat conduction set to " << dom.cont_heat_cond << endl;
      } 
      else 
        cout << "false. "<<endl;      
      } //Rigid bodies
      

    
		
		
		for (auto& ampl : amplitudes) { //TODO: CHECK IF DIFFERENTS ZONES OVERLAP
			// MaterialData* data = new MaterialData();
			int zoneid,valuetype;
			std::vector<double> time, value;
			readValue(ampl["id"], 		zoneid);
			//readValue(zone["valueType"],zoneid);
			readArray(ampl["time"], 	time);
			readValue(ampl["value"], 	value);
			SPH::amplitude amp;
      cout << "Time Size "<<time.size()<<endl;
      amp.id =  zoneid;
			for (int i=0;i<time.size();i++){
				amp.time.push_back(time[i]);
				amp.value.push_back(value[i]);
			}
			dom.amps.push_back(amp);
      //cout << "Amplitude Id "<<zoneid<<" read, "<<i << " values."<<endl;
			//std::cout<< "Zone "<<zoneid<< ", particle count: "<<partcount<<std::	endl;
		}

		for (auto& bc : bcs) { //TODO: CHECK IF DIFFERENTS ZONES OVERLAP
			// MaterialData* data = new MaterialData();
			int zoneid,valuetype,var,ampid;


			double ampfactor = 1.0;
			bool free=true;
			SPH::boundaryCondition bcon;
      bcon.type = 0;        //DEFAULT: VELOCITY
      bcon.valueType = 0;   //DEFAULT: CONSTANT
      bcon.value_ang = 0.0;
			readValue(bc["zoneId"], 	bcon.zoneId);
			readValue(bc["type"], 	bcon.type);
      //type 0 means velocity vc
			readValue(bc["valueType"], 	bcon.valueType);
      readVector(bc["value"], 	      bcon.value);      //Or value linear
      readVector(bc["valueAng"], 	    bcon.value_ang);  //Or Angular value
			if (bcon.valueType == 0){//Constant
        // readVector(bc["value"], 	      bcon.value);      //Or value linear
        // readVector(bc["valueAng"], 	    bcon.value_ang);  //Or Angular value
      } else 
        if ( bcon.valueType == 1){ //Amplitude
				readValue(bc["amplitudeId"], 		bcon.ampId);
				readValue(bc["amplitudeFactor"], 	ampfactor); //DEFAULT NOW IS ZERO
        bcon.ampFactor = ampfactor;
			}
      
      if (bcon.type == Velocity_BC){//Constant
        // readVector(bc["value"], 	      bcon.value);      //Or value linear
        // readVector(bc["valueAng"], 	    bcon.value_ang);  //Or Angular value
      }else if (bcon.type == Convection_BC){
        int count =0;
        double cv_coeff, T_inf;
        readValue(bc["convCoeff"],  bcon.cv_coeff);
        readValue(bc["infTemp"],    bcon.T_inf);
        
        for (size_t a=0; a<dom.Particles.Size(); a++){
          if (dom.Particles[a]->ID == bcon.zoneId) {
            dom.Particles[a]->Thermal_BC	= TH_BC_CONVECTION;
            dom.Particles[a]->T_inf       = bcon.T_inf;
            dom.Particles[a]->h_conv      = bcon.cv_coeff;
            count ++;
          }
        }
        cout << count << " particles have been set with conv coeff: "<<bcon.cv_coeff << "and Tinf "<<bcon.T_inf<<endl;
      } else if (bcon.type == Temperature_BC){
        int count =0;
        readValue(bc["Temp"],    bcon.T);
        //And also applied at every step
        for (size_t a=0; a<dom.Particles.Size(); a++){
          if (dom.Particles[a]->ID == bcon.zoneId) {
            dom.Particles[a]->T	= bcon.T;
            count ++;
          }
        }
        cout << count << " particles have been set with conv coeff: "<<bcon.cv_coeff << "and Tinf "<<bcon.T_inf<<endl;
      }
				
				
			readValue(bc["free"], 	bcon.free);
			dom.bConds.push_back(bcon);
			
      std::cout<< "BCs "<<  ", Zone ID: "<<bcon.zoneId<<", Value :" <<bcon.value<<", Value Type:" <<bcon.valueType<<std::endl;
		}//Boundary Conditions
		
		double IniTemp = 0.;
    int ics_count = 0;
		for (auto& ic : ics){
			double temp;
      int id;
      readValue(ic["Temp"], IniTemp);
      readValue(ic["id"], id);
      cout << "Initial Temp: "<<IniTemp<<endl;

      cout << "Initial condition "<<ics_count<<endl;
			std::vector<double> value(3);
			value[0]=value[1]=value[2]=0.0;
			readValue(ic["value"], 	    value);
      
      int zoneid = -1; //if NO ZONE ID IS ALL SOLID PARTICLES
      readValue(ic["zoneId"], 		zoneid);
      int count = 0, tcount = 0;
      cout << "zone id "<<zoneid<<endl;
			if (dom.Particles.Size()>0)
      if (zoneid == -1){
        cout << "solid part count "<<dom.solid_part_count<<endl;
          cout << "Thermal solver is : "<<dom.thermal_solver<<endl;
        for (size_t a=0; a < dom.solid_part_count; a++){
          dom.Particles[a]->a		= 0.0;
          dom.Particles[a]->v		= Vec3_t(value[0],value[1],value[2]);  

          if (dom.thermal_solver){
            dom.Particles[a]->T = IniTemp;
            tcount ++;
          }
          count ++;
        }
      } else  {
            
          for (size_t a=0; a < dom.Particles.Size(); a++){
            if (dom.Particles[a]->ID == zoneid){
              if (dom.thermal_solver){
                dom.Particles[a]->T = IniTemp;
                count ++;
              }
            }
          }
        }
      
      cout << "Initial condition set wih ID " << id << "applied to "<<count<<" particles."<<endl;
      cout << "TEMP Initial condition set wih ID " << id << "applied to "<<tcount<<" particles."<<endl;
      ics_count++;
		}//ICS
    
    //Add fixed particles, these have priority
    
    //TODO: CHECK IF DIFFERENT ZONES ARE INTERF
    //Generate Domain
    dom.gradKernelCorr = kernel_grad_corr;
    dom.ts_nb_inc = nb_upd_freq;
    
    if (dom.Particles.Size()>0){
    for (size_t a=0; a<dom.Particles.Size(); a++){
      dom.Particles[a]->G				= G;
      dom.Particles[a]->PresEq		= 0;
      dom.Particles[a]->Cs			= Cs;
      dom.Particles[a]->Shepard		= false;
      
      if      ( mattype == "Bilinear" )     dom.Particles[a]->Ep 			= Ep;//HARDENING 
      else if ( mattype == "Hollomon" )     dom.Particles[a]->Material_model  = HOLLOMON;
      else if ( mattype == "JohnsonCook" )  dom.Particles[a]->Material_model  = JOHNSON_COOK;
      else if ( mattype == "GMT" )          dom.Particles[a]->Material_model  = _GMT_;
			dom.Particles[a]->mat             = mat; //NOW MATERIAL IS GIVEN FOR EVERY MATERIAL MODEL (SINCE DAMAGE USES IT)

      dom.Particles[a]->Fail			= 1;
      dom.Particles[a]->Alpha			= alpha;
      dom.Particles[a]->Beta			= beta;
      dom.Particles[a]->TI			= tensins;
      dom.Particles[a]->TIInitDist	= dx;
      dom.Particles[a]->hfac = 1.2; //Only for h update, not used
			
			// if (model_damage)
			// dom.Particles[a]->mat->damage = damage;
      
      // THERMAL PROPS
      dom.Particles[a]->k_T = k_T;
      dom.Particles[a]->cp_T = cp_T;
      dom.Particles[a]->th_exp = th_ex;

			if (mattype == "Hollomon" || mattype == "JohnsonCook" || mattype == "GMT"){ //Link to material is only necessary when it is not bilinear (TODO: change this to every mattype)
       dom.Particles[a]->Sigmay	= mat->CalcYieldStress(0.0,0.0,dom.Particles[a]->T);    
      } else {
        if (Fy>0.0)
          dom.Particles[a]->Sigmay		      = Fy;
        else
          throw new Fatal("Invalid Initial Yield Stress.");
      }

    }
    
    cout << "Reduction Type is: ";
    if (dom.nonlock_sum)
      cout << "Nishimura "<<endl;
    else
      cout << "Locking "<<endl; 
		//dom.SolveDiffUpdateLeapfrog(/*tf*/sim_time,/*dt*/timestep,/*dtOut*/output_time,"test06",1000);
    if (solver=="Mech-Fraser" || solver=="Mech-Thermal-Fraser")
      dom.SolveDiffUpdateFraser(/*tf*/sim_time,/*dt*/timestep,/*dtOut*/output_time,"test06",1e6);
    if (solver=="Mech" || solver=="Mech-Thermal" || solver=="Mech-LeapFrog" || solver=="Mech-Thermal-LeapFrog")
      dom.SolveDiffUpdateLeapFrog(/*tf*/sim_time,/*dt*/timestep,/*dtOut*/output_time,"test06",1e6);
    else if (solver=="Mech-KickDrift" || solver=="Mech-Thermal-KickDrift"){
      dom.SolveDiffUpdateKickDrift(/*tf*/sim_time,/*dt*/timestep,/*dtOut*/output_time,"test06",1e6);
      cout << "Solver-KickDrift is deprecated"<<endl;
    }
    else if (solver=="Thermal")
      dom.ThermalSolve(/*tf*/sim_time,/*dt*/timestep,/*dtOut*/output_time,"test06",1e6);
    else 
      throw new Fatal("Invalid solver.");
		} else {
      throw new Fatal("Particle Count is Null. Please Check Radius and Domain Dimensions.");
    }
		dom.WriteXDMF("maz");
		// dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
		// dom.BC.InOutFlow = 0;

    	//dom.Solve(/*tf*/0.00205,/*dt*/timestep,/*dtOut*/0.00005,"test06",999);
    dom.out_file.close();   
	}	//Argc > 0
  else {cout << "No input file found. Please specify input file."<<endl;}

    return 0;
}

MECHSYS_CATCH
