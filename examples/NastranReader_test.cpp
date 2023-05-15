#include "NastranReader.h"

using namespace SPH;

int main(){
	
	NastranReader reader("Tool.nas");
	reader.WriteCSV("test");
	reader.WriteVTK("test");
	return 0;
	
}