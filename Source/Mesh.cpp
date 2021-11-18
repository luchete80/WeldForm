//#include "Mesh.h"
#include <iostream>
using namespace std;

namespace SPH {

	
TriMesh::TriMesh(){
	
	
}

Element::Element(const int &n1, const int &n2, const int &n3){
	
	//centroid = Vec3_t();
	node[0] = n1; node [1] = n2; node[2] = n3;
	
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
	
	x3 = p1(dir[2]);

	x2=p1(dir[1]); 
	double dl = p(dir[0])/dens;	//Could be allowed 2 diff densities
	//Plane is in 0 and 1 dirs
	for (int j=0; j<dens+1; j++) {
		x1 = p1(dir[0]);
		for (int i=0; i<dens+1; i++){
			Vec3_t v;
			v(dir[0])=x1;v(dir[1])=x2;v(dir[2])=x3;
			node.Push(new Vec3_t(x1,x2,x3));
			cout << "xyz: "<<x1 << ", "<<x2<<", "<<x3<<endl;
			x1+=dl;
		}
		x2+=dl;
	}

	int n[4];
	int el =0;
	
	for (int a = 0 ;a  < 3; a++ ) {
		int j = a;
		for (int i = 0; i < dens; i++ ){
				n[0]=(dens + 1)*j + i; n[1]=n[0] + 1; n[2] = (dens+1)* (j+1) + i ;n[3] = n[2] + 1;
			cout <<" j" << j<<endl;
			int elcon[2][3];	// TODO: check x, y and z normals and node direction 
												// For all plane orientations
			if (positaxisorent) {
				elcon[0][0] = n[0];elcon[0][1] = n[1];elcon[0][2] = n[2];
				elcon[1][0] = n[1];elcon[1][1] = n[3];elcon[1][2] = n[2];
			} else {
				elcon[0][0] = n[0];elcon[0][1] = n[2];elcon[0][2] = n[1];
				elcon[1][0] = n[1];elcon[1][1] = n[2];elcon[1][2] = n[3];				
			}
			cout << "elnodes"<<endl;
			for ( int e= 0; e<2;e++) { // 2 triangles
				element.Push(new Element(elcon[e][0],elcon[e][1],elcon[e][2]));		
				cout << "Element "<< el <<": ";
				for (int en = 0 ; en<3; en++) cout << elcon[e][en]<<", ";
				cout <<endl;
				
				Vec3_t v = ( *node[elcon[e][0]] + *node[elcon[e][1]] + *node[elcon[e][2]] ) / 3. ;
				element[el] -> centroid = v; el++;
				cout << "Centroid" << endl;
			}
		}// i for
		
	}
	///////////////////////////////////////////
	//// MESH GENERATION END

}

//This is done once, Since mesh is rigid
void TriMesh::CalcSpheres(){
	double max;
	cout << "Element radius: "<<endl;
	for (int e = 0; e < element.Size(); e++){ 
		max = 0.;
		Vec3_t rv;
		for (int n = 0 ;n < 3; n++){
			rv = element[e]->node[n] - element[e] -> centroid;
			if (norm(rv) > max) max = norm(rv);
		}
		element[e]->radius = max;
		cout << element[e]->radius<< endl;
	}
	
}

TriMesh::CalcNormals(){
	
}

};