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

#define TAU		0.005
#define VMAX	10.0

#define PRINTVEC(v)	cout << v[0]<<", "<<v[1]<<", "<<v[2]<<endl;

using namespace std;
using namespace SPH;

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
		nlohmann::json bcs 			= j["BoundaryConditions"];
		
		SPH::Domain	dom;
		
		dom.Dimension	= 3;
		
		dom.Nproc	= 4;
		string kernel;
    double ts;
		readValue(config["timeStepSize"], /*scene.timeStepSize*/ts);
    	dom.Kernel_Set(Qubic_Spline);
    
		
		readValue(config["integrationMethod"], dom.Scheme); //0:Verlet, 1:LeapFrog, 2: Modified Verlet

     	//dom.XSPH	= 0.5; //Very important

        double dx,r,h;


		readValue(config["particleRadius"], r);
		double hfactor;
		readValue(config["hFactor"], hfactor);

    	dom.GeneralAfter = & UserAcc;
		
	
		//////////////
		// MATERIAL //
		//////////////
		double rho,E,nu,K,G,Cs,Fy;
    double c[6];
    string mattype;
    cout << "Reading Material.."<<endl;
    readValue(material[0]["type"], 		mattype);
    readValue(material[0]["density0"], 		rho);
    readValue(material[0]["youngsModulus"], 	E);
    readValue(material[0]["poissonsRatio"], 	nu);
    readValue(material[0]["yieldStress0"], 	Fy);
    readValue(material[0]["const1"], 		c[0]);
    readValue(material[0]["const1"], 		c[1]);
    Material_ *mat;
    Elastic_ el(E,nu);
    if (mattype == "Hollomon") mat = new Hollomon(el,Fy,c[0],c[1]);
    //else if (mattype == "JohnsonCook") mat = 
		cout << "Mat type  "<<mattype<<endl;
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
    bool kernel_grad_corr = false;
    readValue(config["artifViscAlpha"],alpha);
    readValue(config["artifViscBeta"],beta);
    readValue(config["contAlgorithm"],cont_alg);
    readValue(config["kernelGradCorr"],kernel_grad_corr);
    dom.auto_ts = auto_ts[0];
    dom.auto_ts_acc = auto_ts[1];
    dom.auto_ts_cont = auto_ts[2];
		////////////
		// DOMAIN //
		////////////
		Vec3_t start,L;
    int id;
		string domtype = "Box";
    int matID;
    bool sym[] = {false,false,false};
		readValue(domblock[0]["id"], 	id);
		readVector(domblock[0]["start"], 	start);
		readVector(domblock[0]["dim"], 	L);
		readValue(domblock[0]["type"], 	domtype); //0: Box
    readValue(domblock[0]["matID"], 	matID); //0: Box
    readBoolVector(domblock[0]["sym"], 	sym); //0: Box
        for (int i=0;i<3;i++) {//TODO: Increment by Start Vector
			dom.DomMax(0) = L[i];
			dom.DomMin(0) = -L[i];
		}		


		
		// inline void Domain::AddCylinderLength(int tag, Vec3_t const & V, double Rxy, double Lz, 
									// double r, double Density, double h, bool Fixed) {
												
		//dom.AddCylinderLength(1, Vec3_t(0.,0.,-L/10.), R, L + 2.*L/10. + dx, r, rho, h, false); 
		
		cout << "Dimensions: "<<endl;
		PRINTVEC(L)
		if (domtype == "Box"){
      cout << "Adding Box Length..."<<endl;      
			dom.AddBoxLength(id ,start, L[0] , L[1],  L[2] , r ,rho, h, 1 , 0 , false, false );		
		}
		else if (domtype == "Cylinder"){
			if (sym[0] && sym[1])
        dom.AddXYSymCylinderLength(0, L[0]/2., L[2], r, rho, h, false, sym[2]); 
      else
        dom.AddCylinderLength(0, start, L[0]/2., L[2], r, rho, h, false, sym[2]); 
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
    if (dim (0)!=0. && dim(1) != 0. && dim(2) !=0. && rigbody_type == "Plane")
      std::cout << "ERROR: Contact Plane Surface should have one null dimension"<<std::endl;
    
    std::vector<TriMesh *> mesh;
    
    if (contact){
      mesh.push_back(new TriMesh);
      mesh[0]->AxisPlaneMesh(2, false, start, Vec3_t(start(0)+dim(0),start(1)+dim(1), dim(2)),40);
      
		}
    
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
			readValue(bc["zoneId"], 	bcon.zoneId);
      //type 0 means velocity vc
			readValue(bc["valueType"], 	bcon.valueType);
			if (bcon.valueType == 0){//Constant
        readVector(bc["value"], 	  bcon.value);
      } else 
        if ( bcon.valueType == 1){ //Amplitude
				readValue(bc["amplitudeId"], 		bcon.ampId);
				readValue(bc["amplitudeFactor"], 	bcon.ampFactor);
			}
				
			readValue(bc["free"], 	bcon.free);
			dom.bConds.push_back(bcon);
			
      std::cout<< "BCs "<<  ", Zone ID: "<<bcon.zoneId<<", Value :" <<bcon.value<<std::endl;
		}
    
    
    //Add fixed particles, these have priority
    
    //TODO: CHECK IF DIFFERENT ZONES ARE INTERF
    //Generate Domain
    dom.gradKernelCorr = kernel_grad_corr;
    if (dom.Particles.Size()>0){
    for (size_t a=0; a<dom.Particles.Size(); a++){
      dom.Particles[a]->G				= G;
      dom.Particles[a]->PresEq		= 0;
      dom.Particles[a]->Cs			= Cs;
      dom.Particles[a]->Shepard		= false;
      
      if ( mattype == "Hollomon" )  dom.Particles[a]->Material_model  = HOLLOMON;
			dom.Particles[a]->mat             = mat;
      dom.Particles[a]->Sigmay		      = Fy;
            
      dom.Particles[a]->Fail			= 1;
      dom.Particles[a]->Alpha			= alpha;
      dom.Particles[a]->Beta			= beta;
      dom.Particles[a]->TI			= 0.3;
      dom.Particles[a]->TIInitDist	= dx;
      dom.Particles[a]->hfac = 1.2;
    }
		dom.SolveDiffUpdateKickDrift(/*tf*/sim_time,/*dt*/timestep,/*dtOut*/output_time,"test06",1000);
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