
// void NastranReader::read( char* fName){
	// string fileName = fName;
	// fstream file;
    // bool found=false;
	// //MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	// cout << "[I] Reading ... "<<endl;
	// file.open(fileName.c_str());
	// if (file.is_open()) {
		// cout << "[I] Found input file " << fileName << endl;
		// found=true;
	// } else {
		// cerr << "[E] Input file " << fileName << " could not be found!!" << endl;
	// }
	
	// int l=0;
	// if (found) {	
		// while(getline(file, line)) {rawData += line + "\n"; l++;}
		// file.close();

		// // Strip all the inline or block comments (C++ style) from the rawData
		// //stripComments(rawData);
		// // Strip all the white spaces from the rawData
		// //strip_white_spaces(rawData);
	// }
	// cout << "[I] "<<l << " lines readed ..." <<endl;

// }