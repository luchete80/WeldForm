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

using namespace std;

class NastranReader {

	std::vector <string> rawData;
	int line_count;
	int element_count;
	
	public:
	NastranReader(char* fName){read(fName);}
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
		cout << "[I] Found input file " << fileName << endl;
		found=true;
	} else {
		cerr << "[E] Input file " << fileName << " could not be found!!" << endl;
	}
	
	int l=0;
	if (found) {	
		while(getline(file, line)) {rawData.push_back(line /*+ "\n"*/); l++;}
		file.close();

		// Strip all the inline or block comments (C++ style) from the rawData
		//stripComments(rawData);
		// Strip all the white spaces from the rawData
		//strip_white_spaces(rawData);
	}
	cout << "[I] "<<l << " lines readed ..." <<endl;
	line_count = l;
	
	//Get number of elements 
	element_count = 0;
	while (!end){
		
	}

}

#endif

