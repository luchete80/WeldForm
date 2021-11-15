#include "Mesh.h"
namespace SPH {

Mesh::Mesh(){
	
	
}

Mesh::AxisPlaneMesh(const int &axis, const int &offset, const Vec3_t p1, const Vec3_t p2,  const int &dens){
	int elemcount = dens * dens;
	
	double x,y,z;
	if (axis == 0 )
	{x = p1(0); y= }
	
	for (int i=0; i<dens+1; i++) {
		for (int j=0; j<dens+1; j++){
			
			node.Push(new Node(x,y,z));
		}
	}

	element.Push(new Element());
	
}

Mesh::CalcNormals(){
	
}

};