//#include "Mesh.h"
namespace SPH {

	
TriMesh::TriMesh(){
	
	
}

Element::Element(const int &n1, const int &n2, const int &n3){
	
	//centroid = Vec3_t();
	
}

// TODO: extend to all dirs
TriMesh::AxisPlaneMesh(const int &axis, bool positaxisorent, const Vec3_t p1, const Vec3_t p2,  const int &dens){
	int elemcount = dens * dens;
	
	double x1,x2,x3;
	double l1,l2;
	Vec3_t p = p2-p1;
	int dir[3];
	if 			(axis == 0 )	{dir[0] = 1; dir[1] = 2;}
	else if (axis == 1 )	{dir[0] = 0; dir[1] = 2;}
	else									{dir[0] = 0; dir[1] = 1;}
	
	dir [2] = axis;
	
	x3 = p(dir[2]);

	x2=p(dir[1]); 
	//Plane is in 0 and 1 dirs
	for (int j=0; j<dens+1; j++) {
		x1 = p(dir[0]);
		for (int i=0; i<dens+1; i++){
			Vec3_t v;
			v(dir[0])=x1;v(dir[1])=x2;v(dir[2])=x3;
			node.Push(new Vec3_t(x1,x2,x3));
			x1+=dens;
		}
		x2+=dens;
	}

	int n[4];
	int el =0;

	for (int j = 0 ; j < dens; j++ ) {
		for (int i=0; i < dens; i++ ){
				n[0]=i; n[1]=i+1; n[2] = dens * j;n[3] = dens * j + 1;
			
			int elcon[2][3];	// TODO: check x, y and z normals and node direction 
												// For all plane orientations
			if (positaxisorent) {
				elcon[0][0] = n[0];elcon[0][1] = n[1];elcon[0][2] = n[2];
				elcon[1][0] = n[1];elcon[1][1] = n[3];elcon[1][2] = n[2];
			} else {
				elcon[0][0] = n[0];elcon[0][1] = n[2];elcon[0][2] = n[1];
				elcon[1][0] = n[1];elcon[1][1] = n[2];elcon[1][2] = n[3];				
			}
			
			for ( int e= 0; e<2;e++) { // 2 triangles
				element.Push(new Element(elcon[e][0],elcon[e][1],elcon[e][2]));			
				Vec3_t v = ( *node[elcon[e][0]] + *node[elcon[e][1]] + *node[elcon[e][2]] ) / 3. ;
				element[el] -> centroid = v; el++;
			}

		}
	}
	///////////////////////////////////////////
}

TriMesh::CalcNormals(){
	
}

};