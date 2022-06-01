#ifndef _PLANE_H_
#define _PLANE_H_

#include "matvec.h"

namespace SPH{
 
class Plane{
public:

  	double 	pplane;			//In boundary elements, plane coefficient, useful for contact
    Vec3_t 	normal;
};

};


#endif