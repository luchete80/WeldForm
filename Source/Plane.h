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
  	double 	pplane;			//In boundary elements, plane coefficient, useful for contact
    Vec3_t 	normal;
    Vec3_t 	tg[2];
};

};


#endif