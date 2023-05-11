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
				if (domi.bConds[bc].type == 0 ){ //VELOCITY
					domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);
					domi.Particles[i]->v		= domi.bConds[bc].value;
				}
			}
			
		}
	}
  
  if (domi.contact){
    for (int bc=0;bc<domi.bConds.size();bc++){
      for (int m=0;m<domi.trimesh.size();m++){
        if (domi.trimesh[m]->id == domi.bConds[bc].zoneId)
          domi.trimesh[m]->SetVel(domi.bConds[bc].value);
      }//mesh
    }//bcs
  }//contact
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
		
		dom.Dimension	= 3;
		
		string kernel;
    double ts;
    
    cout << "--------------------------------------------"<<endl;
    cout << "----------------- WELDFORM -----------------"<<endl;
    cout << "----------------- v. 0.4.1 -----------------"<<endl;
    cout << "--------------------------------------------"<<endl<<endl<<endl;
    
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
    
    if (solver=="Mech-Thermal")
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
    double Et, Ep;  //Hardening (only for bilinear and multilear)
    std::vector<double> c;
    c.resize(10);
    string mattype = "Bilinear";
    cout << "Reading Material.."<<endl;
    cout << "Type.."<< endl; readValue(material[0]["type"], 		mattype);
    cout << "Density.."<< endl; readValue(material[0]["density0"], 		rho);
    readValue(material[0]["youngsModulus"], 	E);
    readValue(material[0]["poissonsRatio"], 	nu);
    readValue(material[0]["yieldStress0"], 	Fy);
    readArray(material[0]["const"], 		c);
    Material_ *mat;
    Elastic_ el(E,nu);
    cout << "Mat type  "<<mattype<<endl;
    if      (mattype == "Bilinear")    {
      Ep = E*c[0]/(E-c[0]);		                              //only constant is tangent modulus
      cout << "Material Constants, Et: "<<c[0]<<endl;
    } else if (mattype == "Hollomon")    {
      mat = new Hollomon(el,Fy,c[0],c[1]);
      cout << "Material Constants, K: "<<c[0]<<", n: "<<c[1]<<endl;
    } else if (mattype == "JohnsonCook") {
      //Order is 
                                 //A(sy0) ,B,  ,C,   m   ,n   ,eps_0,T_m, T_transition
      mat = new JohnsonCook(el,Fy, c[0],c[1],c[3],c[2],c[6], c[4],c[5]); //First is hardening // A,B,C,m,n_,eps_0,T_m, T_t);	 //FIRST IS n_ than m
      cout << "Material Constants, B: "<<c[0]<<", C: "<<c[1]<<", n: "<<c[2]<<", m: "<<c[3]<<", T_m: "<<c[4]<<", T_t: "<<c[5]<<", eps_0: "<<c[6]<<endl;
    } else                              throw new Fatal("Invalid material type.");
    
    
    // THERMAL PROPERTIES

    double k_T, cp_T;
    readValue(material[0]["thermalCond"], 	  k_T);
    readValue(material[0]["thermalHeatCap"], 	cp_T);    
    
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
    string cont_alg = "Fraser";
    bool auto_ts[] = {true, false, false}; //ONLY VEL CRITERIA
		readValue(config["cflMethod"], cflMethod);
		if (cflMethod == 0)
			readValue(config["timeStepSize"], timestep);
		else {
			readValue(config["cflFactor"], cflFactor);
			timestep = (cflFactor*h/(Cs));
		}
    readValue(config["outTime"], output_time);
    readValue(config["simTime"], sim_time);
    readBoolVector(config["autoTS"], auto_ts);
    double alpha = 1.;
    double beta = 0.;
    bool h_upd = false;
    double tensins = 0.3;
    bool kernel_grad_corr = false;
    readValue(config["artifViscAlpha"],alpha);
    readValue(config["artifViscBeta"],beta);
    readValue(config["contAlgorithm"],cont_alg);
    readValue(config["kernelGradCorr"],kernel_grad_corr);
    readValue(config["smoothlenUpdate"],h_upd);
    dom.auto_ts = auto_ts[0];
    dom.auto_ts_acc = auto_ts[1];
    dom.auto_ts_cont = auto_ts[2];
		
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
    bool sym[] = {false,false,false};
		readValue(domblock[0]["id"], 	id);
		readVector(domblock[0]["start"], 	start);
		cout << "Reading Domain dim" << endl;  readVector(domblock[0]["dim"], 	L);
		cout << "Reading Domain type" << endl; readValue(domblock[0]["type"], 	domtype); //0: Box
    cout << "Reading Domain mat id" << endl;  readValue(domblock[0]["matID"], 	matID); //0: Box
    cout << "Grid Coordinate System" << endl;  readValue(domblock[0]["gridCoordSys"], 	gridCS); //0: Box
    readBoolVector(domblock[0]["sym"], 	sym); //0: Box
        for (int i=0;i<3;i++) {//TODO: Increment by Start Vector
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
			dom.AddBoxLength(id ,start, L[0] , L[1],  L[2] , r ,rho, h, 1 , 0 , false, false );		
		}
		else if (domtype == "Cylinder"){
      cout << "Adding Cylinder";      
			if (sym[0] && sym[1]){
        cout << " with symmetry..."<<endl;
        dom.AddXYSymCylinderLength(0, L[0]/2., L[2], r, rho, h, false, sym[2]); 
      }
      else {
        cout << "..."<<endl;
        if ( gridCS == "Cartesian")
          dom.AddCylinderLength(0, start, L[0]/2., L[2], r, rho, h, false, sym[2]); 
        else if (gridCS == "Cylindrical")
          dom.AddCylUniformLength(0, L[0]/2.,L[2], r, rho, h);
          
      }
    }

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
    string rigbody_type;
    bool contact = false;
    if (readValue(rigbodies[0]["type"],rigbody_type))
      contact = true;
    Vec3_t dim;
    
		readVector(rigbodies[0]["start"], 	start);       
		readVector(rigbodies[0]["dim"], 	dim); 
    bool flipnormals = false;
    readValue(rigbodies[0]["flipNormals"],flipnormals);
    
    double heatcap = 1.;
    readValue(rigbodies[0]["thermalHeatCap"],heatcap);
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
    std::vector<TriMesh *> mesh;
    
    cout << "Set contact to ";
    if (contact){
      cout << "true."<<endl;
      dom.contact = true;
      cout << "Reading contact mesh..."<<endl;
      //TODO: CHANGE TO EVERY DIRECTION
      int dens = 10;
      readValue(rigbodies[0]["partSide"],dens);
      if (rigbody_type == "Plane"){
        // TODO: CHECK IF MESH IS NOT DEFINED
        mesh.push_back(new TriMesh);
        mesh[0]->AxisPlaneMesh(2, false, start, Vec3_t(start(0)+dim(0),start(1)+dim(1), start(2)),dens);
      } else if (rigbody_type == "File"){
        string filename = "";
        readValue(rigbodies[0]["fileName"], 	filename); 
        cout << "Reading Mesh input file " << filename <<endl;
        SPH::NastranReader reader(filename.c_str());
          mesh.push_back (new SPH::TriMesh(reader,flipnormals ));
      }
      cout << "Creating Spheres.."<<endl;
      //mesh.v = Vec3_t(0.,0.,);
      mesh[0]->CalcSpheres(); //DONE ONCE
      double hfac = 1.1;	//Used only for Neighbour search radius cutoff
      cout << "Adding mesh particles ...";
      int id;
      readValue(rigbodies[0]["zoneId"],id);
      dom.AddTrimeshParticles(mesh[0], hfac, id); //AddTrimeshParticles(const TriMesh &mesh, hfac, const int &id){
        
      
      std::vector<double> fric_sta(1), fric_dyn(1), heat_cond(1);
      readValue(contact_[0]["fricCoeffStatic"], 	fric_sta[0]); 
      readValue(contact_[0]["fricCoeffDynamic"], 	fric_dyn[0]); 
      readValue(contact_[0]["heatCondCoeff"], 	  heat_cond[0]);
      
      bool heat_cond_ = false;
      if (readValue(contact_[0]["heatConductance"], 	heat_cond_)){
        dom.cont_heat_cond = true;
        dom.contact_hc = heat_cond[0];
      }
      
      dom.friction_dyn = fric_dyn[0];
      dom.friction_sta = fric_sta[0];
      dom.PFAC = 0.8;
      dom.DFAC = 0.0;
      
		} 
    else 
      cout << "false. "<<endl;
    
		std::vector <SPH::amplitude> amps;
		
		for (auto& ampl : amplitudes) { //TODO: CHECK IF DIFFERENTS ZONES OVERLAP
			// MaterialData* data = new MaterialData();
			int zoneid,valuetype;
			std::vector<double> time, value;
			readValue(ampl["id"], 		zoneid);
			//readValue(zone["valueType"],zoneid);
			readArray(ampl["time"], 	time);
			readValue(ampl["value"], 	value);
			SPH::amplitude amp;
			for (int i=0;i<time.size();i++){
				amp.time.push_back(time[i]);
				amp.value.push_back(value[i]);
			}
			amps.push_back(amp);
			//std::cout<< "Zone "<<zoneid<< ", particle count: "<<partcount<<std::	endl;
		}

		for (auto& bc : bcs) { //TODO: CHECK IF DIFFERENTS ZONES OVERLAP
			// MaterialData* data = new MaterialData();
			int zoneid,valuetype,var,ampid;

			double ampfactor;
			bool free=true;
			SPH::boundaryCondition bcon;
      bcon.type = 0;        //DEFAULT: VELOCITY
      bcon.valueType = 0;   //DEFAULT: CONSTANT
      bcon.value_ang = 0.0;
			readValue(bc["zoneId"], 	bcon.zoneId);
      //type 0 means velocity vc
			readValue(bc["valueType"], 	bcon.valueType);
			if (bcon.valueType == 0){//Constant
        readVector(bc["value"], 	      bcon.value);      //Or value linear
        readVector(bc["valueAng"], 	    bcon.value_ang);  //Or Angular value
      } else 
        if ( bcon.valueType == 1){ //Amplitude
				readValue(bc["amplitudeId"], 		bcon.ampId);
				readValue(bc["amplitudeFactor"], 	bcon.ampFactor);
			}
				
			readValue(bc["free"], 	bcon.free);
			dom.bConds.push_back(bcon);
			
      std::cout<< "BCs "<<  ", Zone ID: "<<bcon.zoneId<<", Value :" <<bcon.value<<std::endl;
		}//Boundary Conditions
		
		double IniTemp = 0.;
		for (auto& ic : ics){
			double temp;
			if (solver == "Mech-Thermal"){
				readValue(ic["Temp"], IniTemp);
				cout << "Initial Temp: "<<IniTemp<<endl;
			}
		}
    
    //Add fixed particles, these have priority
    
    //TODO: CHECK IF DIFFERENT ZONES ARE INTERF
    //Generate Domain
    dom.gradKernelCorr = kernel_grad_corr;
    dom.ts_nb_inc = 5;
    
    if (dom.Particles.Size()>0){
    for (size_t a=0; a<dom.Particles.Size(); a++){
      dom.Particles[a]->G				= G;
      dom.Particles[a]->PresEq		= 0;
      dom.Particles[a]->Cs			= Cs;
      dom.Particles[a]->Shepard		= false;
      
      if      ( mattype == "Bilinear" )     dom.Particles[a]->Ep 			= Ep;//HARDENING 
      else if ( mattype == "Hollomon" )     dom.Particles[a]->Material_model  = HOLLOMON;
      else if ( mattype == "JohnsonCook" )  dom.Particles[a]->Material_model  = JOHNSON_COOK;
			if (mattype == "Hollomon" || mattype == "JohnsonCook"){ //Link to material is only necessary when it is not bilinear (TODO: change this to every mattype)
        dom.Particles[a]->mat             = mat;
       dom.Particles[a]->Sigmay	= mat->CalcYieldStress(0.0,0.0,0.0);    
      }
      dom.Particles[a]->Sigmay		      = Fy;
            
      dom.Particles[a]->Fail			= 1;
      dom.Particles[a]->Alpha			= alpha;
      dom.Particles[a]->Beta			= beta;
      dom.Particles[a]->TI			= tensins;
      dom.Particles[a]->TIInitDist	= dx;
      dom.Particles[a]->hfac = 1.2; //Only for h update, not used
      
      // THERMAL PROPS
      dom.Particles[a]->k_T = k_T;
      dom.Particles[a]->cp_T = cp_T;
	  
	  dom.Particles[a]->T = IniTemp;
    }
    
    cout << "Reduction Type is: ";
    if (dom.nonlock_sum)
      cout << "Nishimura "<<endl;
    else
      cout << "Locking "<<endl;
		//dom.SolveDiffUpdateLeapfrog(/*tf*/sim_time,/*dt*/timestep,/*dtOut*/output_time,"test06",1000);
    if (solver=="Mech" || solver=="Mech-Thermal")
      dom.SolveDiffUpdateFraser(/*tf*/sim_time,/*dt*/timestep,/*dtOut*/output_time,"test06",1000);
    else if (solver=="Mech" || solver=="Mech-Thermal-KickDrift")
      dom.SolveDiffUpdateKickDrift(/*tf*/sim_time,/*dt*/timestep,/*dtOut*/output_time,"test06",1000);
    else if (solver=="Thermal")
      dom.ThermalSolve(/*tf*/sim_time,/*dt*/timestep,/*dtOut*/output_time,"test06",1000);
    else 
      throw new Fatal("Invalid solver.");
		} else {
      throw new Fatal("Particle Count is Null. Please Check Radius and Domain Dimensions.");
    }
		dom.WriteXDMF("maz");
		// dom.m_kernel = SPH::iKernel(dom.Dimension,h);	
		// dom.BC.InOutFlow = 0;

    	//dom.Solve(/*tf*/0.00205,/*dt*/timestep,/*dtOut*/0.00005,"test06",999);
	}	//Argc > 0
  else {cout << "No input file found. Please specify input file."<<endl;}
	
    return 0;
}

MECHSYS_CATCH
