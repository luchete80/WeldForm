#include "matvec.h" 

namespace SPH {
	
void Domain::AddTrimeshParticles(const TriMesh &mesh, const int &id){
	
	double Density =0.;
	double h = 1.1;
	bool Fixed = false;
	 
	for ( int e = 0; e < mesh.element.Size(); e++ ){
		Vec3_t pos = mesh.element[e]->centroid;
		Particles.Push(new Particle(id,pos,Vec3_t(0,0,0),0.0,Density,h,Fixed));
	}
}


};