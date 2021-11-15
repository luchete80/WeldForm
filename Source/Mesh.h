namespace SPH{
// It is like triangular mesh
class Sphere{
	public:
	Vec3_t pos;
	double radius;
};

class Element {
	int node[3];			//Node position (clockwise outer normal)
	Sphere centroid;	//With radius equal to size
	
};
class Mesh{
	
	public:

	Array <Vec3_t*> node;
	Array <Element*> element;
	Mesh();
	AxisPlaneMesh(const int &axis, const int &offset, const Vec3_t p1, const Vec3_t p2, const int &dens);
	CalcNormals();
};

};