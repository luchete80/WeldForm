#ifndef NASTRAN_READER_H_
#define NASTRAN_READER_H_

#include <map>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <cctype>
#include <string>
#include <iostream>
#include <fstream>

#include "Mesh.h"

using namespace std;

class NastranReader {

	std::vector <string> rawData;
	int line_count;
	int element_count;
  
  //Flattened arrays such as GPU type in order of mantain this
  double  *node;
  int     *elcon;
  
  //TriMesh trimesh;
	
	public:
	NastranReader(char* fName){read(fName);}
  ~NastranReader();
	inline void read(char *fName);
	
};

void NastranReader::read( char* fName){
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
  int node_count = 0;
  int elem_count = 0;
  
  bool start_node = false;
  bool start_elem = false;
  
  int line_start_node;
  
	if (found) {	
		while(getline(file, line)) {
      rawData.push_back(line /*+ "\n"*/); 
      //if (strcmp(str_inp1, str_inp2) == 0)
      //Increment nodes
      //if (strcmp(str_inp1, str_inp2) == 0)
        //or str.compare
      //cout << "Searching "<<line.substr(0,4)<<endl;
      if (line.substr(0,4) == string("GRID")){
        //cout << "Node found!"<<endl;
        node_count++;
        if (!start_node){
          start_node = true;
          line_start_node = l;
        }
      } else if (line.substr(0,5) == string("CTRIA")){
        //cout << "Element found!"<<endl;
        elem_count++;
      }
      l++;
    }
		file.close();
    
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
  node  = new double [3 * node_count];
  int curr_line = line_start_node;
  for (int n=0;n<node_count;n++){
    string temp = rawData[curr_line].substr(0,5)
  }
  
  //IF FIXED FIELD
  cout << "Allocating Elements..."<<endl;
  elcon = new int    [3 * elem_count];
    
  cout << "Done."<<endl;
	
}

NastranReader::~NastranReader(){
  
  delete node;
  delete elcon;  
}

#endif

