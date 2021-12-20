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
	Vec3_t 	normal;
	Vec3_t 	v;					//At centroid
	double 	radius;
	int 		node[3];
	double 	pplane;			//In boundary elements, plane coefficient, useful for contact
	int 		nfar;						//farthest away node from baricenter
	//Sphere* centroid;
	//Mesh*		mesh;
};


class TriMesh{
	
	public:

	Array <Element* > 	element;
	Array <Vec3_t* > 		node;
	Array <Vec3_t* > 		node_v;				//Node velocities
	
	Vec3_t							v;						//Constant Uniform v
	TriMesh();
	inline void AxisPlaneMesh(const int &axis, bool positaxisorent, const Vec3_t p1, const Vec3_t p2, const int &dens);
	inline void ApplyConstVel(const Vec3_t &v);
	inline void CalcCentroidVelFromNodes();
	inline void UpdatePlaneCoeff();
	inline void UpdatePos(const double &dt);
	inline void CalcNormals();
	inline void CalcSpheres();
	void CalcCentroids();
};

};
#include "Mesh.cpp"

#endif