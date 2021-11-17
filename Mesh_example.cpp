#include "Mesh.h"
#include <iostream>

using namespace SPH;
using namespace std;

int main(){
	
	TriMesh mesh;
	
	cout << "Creating Mesh" << endl;
	mesh.AxisPlaneMesh(2,true,Vec3_t(0.,0.,0.),Vec3_t(1.,1.,0.),10);
	
	
}