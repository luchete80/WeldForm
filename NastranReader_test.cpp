#include "NastranReader.h"



int main(){
	
	NastranReader reader("Tool.nas");
	reader.WriteCSV("test");
	reader.WriteVTK("test");
	return 0;
	
}