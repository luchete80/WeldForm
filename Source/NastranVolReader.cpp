//#include "NastranVolReader.h"
#include "Plane.h"

namespace SPH{
  
void NASElement::CalcCentroids() {
    centroid = Vec3_t(0.,0.,0.);
    for (int i=0;i<nodecount;i++) centroid += *mesh->nod[node[i]];
    centroid /= nodecount;
  }

class Tri {
  public:
  
  Tri(Vec3_t &v1,Vec3_t &v2,Vec3_t &v3){
    v[0] = v1; v[1] = v2; v[2] = v3;
    for (int i=0;i<3;i++){
      int k = i + 1;
      if (k == 3) k = 0; 
      l[i] = norm(v[k]-v[i]);
    }
    p = 0.5*(l[0]+l[1]+l[2]);
    area = sqrt(p*(p-l[0])*(p-l[1])*(p-l[2]));
  }
  Vec3_t v[3];
  double l[3];
  double area;
  double p;
};

void Tetra::CalcVol() {
    //1/3 AB .h 
    
    //AB = sqrt (p (p-a)(p-b)(p-c)) //HERON formula Where p is (a+b+c)/2
    Tri tri(*mesh->nod[node[0]],*mesh->nod[node[1]],*mesh->nod[node[2]]);

    Plane plane(*mesh->nod[node[0]],*mesh->nod[node[1]],*mesh->nod[node[2]]);
    double h = plane.dist_to_point(*mesh->nod[node[3]]);
    
    
    //cout << "test X+h.n: "<< dot(*mesh->nod[node[3]] - h*normal, normal) << ", d "<<d<<endl; 
    vol = 1./3. * tri.area * h;
}

//Nodes are CCWW from bottom to top
void Hexa::CalcVol(){
  Plane plane(*mesh->nod[node[0]],*mesh->nod[node[1]],*mesh->nod[node[2]]);  
  Tri tri(*mesh->nod[node[0]],*mesh->nod[node[1]],*mesh->nod[node[2]]);
  Tri tri2(*mesh->nod[node[1]],*mesh->nod[node[2]],*mesh->nod[node[3]]);
  double ab = tri.area + tri2.area;
  //Get midpoint of top plane 
  
  //Projecting top verices to plane 1
  
  Vec3_t vmed = 1./4.*(*mesh->nod[node[4]]+*mesh->nod[node[5]]+*mesh->nod[node[6]]+*mesh->nod[node[7]]);
  double h = plane.dist_to_point(vmed);
  double vol1 = ab * h;
  Plane plane_top(plane.normal, vmed);
  Vec3_t vtop[4];
  for (int i=0;i<4;i++){
    vtop[i] = plane_top.project_point(*mesh->nod[node[4+i]]);
  }
  
  Tri tri3(vtop[0],vtop[1],vtop[2]);
  Tri tri4(vtop[1],vtop[2],vtop[3]);
  
  double ab2 = tri3.area + tri4.area;
  
  vol = 0.5*(ab+ab2)*h;
  //cout << "Vol "<<vol <<endl; 
  // TODO: Do the same as opposite with top plane and average! 
  
  
  //cout << "d, X.n "<<plane.pplane<< ", " << dot(plane.normal, *mesh->nod[node[0]])<< endl;
  //cout << "normals "<<plane.normal<< ", " <<plane2.normal<<endl;
}

//Always  seems to be nodes 3 and 6 to be collapsing
//TODO: EVALUATE ALL CASES
//Order are CCW and from bottom to top

void Pyram::CalcVol(){
  //So square is at nodes 0,1, 3 and 4
  // 3 4
  // 0 1
  Tri tri(*mesh->nod[node[0]],*mesh->nod[node[1]],*mesh->nod[node[3]]);
  Tri tri2(*mesh->nod[node[0]],*mesh->nod[node[1]],*mesh->nod[node[4]]);  

  Plane plane(*mesh->nod[node[0]],*mesh->nod[node[1]],*mesh->nod[node[3]]);
    
  double h = plane.dist_to_point(*mesh->nod[node[2]]);
  
  
  //cout << "test X+h.n: "<< dot(*mesh->nod[node[3]] - h*normal, normal) << ", d "<<d<<endl; 
  vol = 1./3. * (tri.area + tri2.area) * h;
  //cout << "vol "<<vol<<endl;
  
}
  
void NastranVolReader::read( const char* fName){
	string fileName = fName;
  string line;
  //rawData="";
	fstream file;
    bool found=false;
	//MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	cout << "[I] Reading ... "<<endl;
	file.open(fileName.c_str());
	if (file.is_open()) {
		//cout << "[I] Found input file " << fileName << endl;
		found=true;
	} else {
		//cerr << "[E] Input file " << fileName << " could not be found!!" << endl;
	}
	
	int l=0;
  node_count = 0;
  elem_count = 0;
  
  bool start_node = false;
  bool start_elem = false;
  
  int line_start_node;
	int line_start_elem;
  bool continuation_line; //Assumes only 2 lines
	if (found) {	
		while(getline(file, line)) {
      //if (strcmp(str_inp1, str_inp2) == 0)
      //Increment nodes
      //if (strcmp(str_inp1, str_inp2) == 0)
        //or str.compare
      //cout << "Searching "<<line.substr(0,4)<<endl;
      if (line.substr(0,4) == string("GRID")){
        rawData.push_back(line /*+ "\n"*/); 
        //cout << "Node found!"<<endl;
        node_count++;
        if (!start_node){
          start_node = true;
          line_start_node = l;
        }
        l++;
      } else if (line.substr(0,6) == string("CTETRA")){
        rawData.push_back(line /*+ "\n"*/); 
        if (!start_elem){
          start_elem = true;
					line_start_elem = l;
          cout << "Tetra found"<<endl;
        }
        //cout << "Element found!"<<endl;
        //cout << "CTETRA"<<endl;
        elem_count++;
        l++;
      } else if (line.substr(0,5) == string("CHEXA")){
        continuation_line = true;
        rawData.push_back(line /*+ "\n"*/); 
        if (!start_elem){
          start_elem = true;
					line_start_elem = l;
          cout << "Hexa found at line "<< l << endl;
				}
        elem_count++;
        l++;
      } else if (line.substr(0,6) == string("CPENTA")){
        rawData.push_back(line /*+ "\n"*/); 
        if (!start_elem){
          start_elem = true;
					line_start_elem = l;
          cout << "Penta found"<<endl;
				}
        elem_count++;
        l++;
      } else if (line.substr(0,1) == string("+")){
        rawData.push_back(line /*+ "\n"*/);         
      }
      else if (continuation_line){ // comments or continuation 
        
      }
      //cout << "Line "<<l<<endl;
    }
		file.close();
    cout << "Line start elem "<<line_start_elem<<endl;
    
    cout << node_count <<" nodes and "<<elem_count<< " elements found."<<endl;

		// Strip all the inline or block comments (C++ style) from the rawData
		//stripComments(rawData);
		// Strip all the white spaces from the rawData
		//strip_white_spaces(rawData);
	}
	cout << "[I] "<<l << " lines readed ..." <<endl;
	line_count = l;
  
  //Allocating nodes 
  cout << "Allocating nodes"<<endl;
  //node  	= new double 	[3 * node_count];
  

	// NODAL FIELD DATA IS: GRID|ID|CP|X1|	
  int curr_line = line_start_node;
	l = curr_line;
	Vec3_t min( 1000., 1000., 1000.);
  Vec3_t max(-1000.,-1000.,-1000.);
	
  std::vector <int> nodeid(node_count);
  
	for (int n=0;n<node_count;n++){
    //cout << n+1; //DEBUG
		string temp = rawData[l].substr(FIELD_LENGTH,FIELD_LENGTH); //Second field, id
		nodeid[n] = atoi(temp.c_str());
		nodepos.insert(std::make_pair(atoi(temp.c_str()),n));
		//cout << "id: "<<nodeid[n]<<endl;
    double ncoord[3];
		for (int i = 0;i<3;i++) {
			int pos = 3*(FIELD_LENGTH)+ i*FIELD_LENGTH;
			//cout << "pos: "<<pos<<endl; 
			string temp = rawData[l].substr(pos,FIELD_LENGTH);
			
			//Search + or - symbol (SCIENTIFIC NOTATION)
			//BUT! If these signs (mainly the minus), are the FIRST character there is 
			//not the exponential sign
			int sign_pos = 0;
			if 			(temp.find("+")!=std::string::npos) {
				sign_pos = temp.find("+"); //Performance....
				temp.insert(sign_pos,"E");
			}
			else if (temp.find("-")!=std::string::npos) {
				sign_pos = temp.find("-");
				if (sign_pos!=0)
					temp.insert(sign_pos,"E");
				else { //Search for another "-" (NOW IT WILL BE SCIENTIFIC NOTATION)
					sign_pos = temp.find("-",1);
					if (sign_pos !=std::string::npos) {
						temp.insert(sign_pos,"E");
					}
				}
			}	
			
			double d = strtod(temp.c_str(),NULL);
      ncoord[i] = d;
			//cout << temp<<", conv: "<<d<<"sign pos" << sign_pos<<endl;
			//cout <<d<< " ";
			//node[3*n+i] = d;
			if (d<min[i])
				min[i] = d;
			else if (d > max[i])
				max[i] = d;
		}
    nod.push_back(new Vec3_t(ncoord[0],ncoord[1],ncoord[2]));
		l++;
  }
	
	cout << "Min values: "<< min <<endl;
	cout << "Max values: "<< max <<endl;	
  
  //IF FIXED FIELD
  cout << "Allocating Elements..."<<endl;
	// ASSUMING NODE IS FROM 1
  //elcon = new int    [3 * elem_count];
  int cypracount =0;
  int chexacount, ctetracount, cprismcount;  
  chexacount = ctetracount = cprismcount = 0;
	map<int, int>::iterator it;
  curr_line = line_start_elem;
	l = curr_line;
  int line_incr;
  for (int n=0;n<elem_count;n++){
    std::vector<int> nv;
    int nodecount;
    line_incr = 1;
    //cout << n+1<< " ";
    if (rawData[l].substr(0,6) == string("CTETRA")){
      //cout << "penta "<<endl;
      nodecount = 4;
      ctetracount++;
    }
    if (rawData[l].substr(0,6) == string("CPENTA")){ //All in one ROW if prism, what about pyramid?
      //cout << "penta "<<endl;
      nodecount = 6;
    }
    else if (rawData[l].substr(0,5) == string("CHEXA")){
      //cout << "penta "<<endl;
      line_incr = 2;
      nodecount = 8;
      chexacount++;
    }
    // if (rawData[l].substr(0,6) == string("CPYRAM")){ //All in one ROW if prism, what about pyramid?
      // //cout << "penta "<<endl;
      // nodecount = 5;
    // }

    int linecount = 0;
    int nid[8];
    bool ispyra = false;
		for (int en=0;en<nodecount;en++){
      int lnum = 0;
      string temp;
      int pos, d;
      if ( en < 6 ) {//MAX FIELDS PER ROW 
        pos = 3*(FIELD_LENGTH)+ en*FIELD_LENGTH; //First 3 are ELEMTYE, el id  and PROP
        nid[en]=pos;
			} else {
        pos = 1*(FIELD_LENGTH)+ (en-6)*FIELD_LENGTH; //First 3 are ELEMTYE, el id  and PROP
        lnum = 1; 
      }
			temp = rawData[l+lnum].substr(pos,FIELD_LENGTH); //Second field, id
      //cout << "Field "<<temp<<endl;
      d = atoi(temp.c_str());
      nid[en] = d;
      int nod = nodepos.find(d)->second;
      
      if (nodecount == 6) {//Check for collapsing nodes (CPYARM)
        for (int n=0;n<en;n++){
          if (d == nid[n]){
            cypracount++;
            ispyra = true;
          }
        }
      }
      if (!ispyra) nv.push_back(nod);
     // cout <<d<<", "<<nod<<"; "; 
      
			//elcon[3*n+en] = nod;
			//cout << d<<" ";
		}//for node count;
    if (!ispyra &&nodecount ==6 )cprismcount++;
    //cout << endl;
    if (nv.size() == 4)
      elem.push_back(new Tetra(this,nv));
    else if (nv.size() == 5)
      elem.push_back(new Pyram(this,nv));
    else if (nv.size() == 6)
      elem.push_back(new Prism(this,nv));
    else if (nv.size() == 8)
      elem.push_back(new Hexa(this,nv));
		//cout << endl;
		l += line_incr;
	}// For elem count    
  cout << "Done."<<endl;
  cout << "Element count "<<endl;
  cout << "Tetra: "<<ctetracount<<endl;
  cout << "Hexa: "<<chexacount<<endl;
  cout << "Pyramid: "<<cypracount<<endl;
  cout << "Prism: "<<cprismcount<<endl;
	
}

}; //SPH