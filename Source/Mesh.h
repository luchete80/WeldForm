#ifndef _MESH_H_
#define _MESH_H_

#include "matvec.h"
//#include

namespace SPH{
	
// class Sphere{
	// public:
	// Sphere();
	// Sphere (Element *elem);
	
	// Vec3_t  pos;
	// //Vec3_t* x[3];
	// Element *element;
	// int node[3];	//CONST 
	// double radius;
	
// };


// It is like triangular mesh
class Element{
	public:
	Element(){}
	Element(const int &n1, const int &n2, const int &n3);
	
	//SPHERE
	Vec3_t 	centroid;
	double 	radius;
	int 		node[3];
	//Sphere* centroid;
	//Mesh*		mesh;
};


class TriMesh{
	
	public:

	Array <Element* > element;
	Array <Vec3_t* > node;
	TriMesh();
	inline void AxisPlaneMesh(const int &axis, bool positaxisorent, const Vec3_t p1, const Vec3_t p2, const int &dens);
	inline void CalcNormals();
	inline void CalcSpheres();
};

};
#include "Mesh.cpp"

#endif