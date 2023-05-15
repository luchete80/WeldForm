#ifndef _PLANE_H_
#define _PLANE_H_

#include "matvec.h"

namespace SPH{
 
class Plane{
public:
    Plane(){}
    Plane (Vec3_t &n, Vec3_t &tg0, Vec3_t &tg1, double &pp){
      normal=n;
      pplane=pp;
      tg[0] = tg0;
      tg[1] = tg1;
    }
    Plane (Vec3_t &p1,Vec3_t &p2,Vec3_t &p3){
    	Vec3_t w;
      
      p = p1;
      tg[0] = p3 - p1;
      tg[1] = p2 - p1;
      w = cross(tg[0],tg[1]);
      //cout << "u, v, w "<<u<<", " <<v << ", " <<w<<endl;
      normal = w/norm(w);
      pplane = dot (p1,normal);      
    }
    Plane (Vec3_t &n, Vec3_t &p_){
      //Normalize!
      normal = n/norm(n);
      p = p_;
      pplane = dot(n,p);
      
    }
    double dist_to_point(Vec3_t q){
      return abs(dot (q-p,normal));
    }
    Vec3_t project_point(Vec3_t &q){
      double h = dist_to_point(q);
      Vec3_t r = q - h * normal;
      return r;
    }
    
  	double 	pplane;			//In boundary elements, plane coefficient, useful for contact
    Vec3_t 	normal;
    Vec3_t 	tg[2];
    Vec3_t p;      //Necesary to distance to plane
};

};


#endif