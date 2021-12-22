#ifndef _MATERIAL_H_
#define _MATERIAL_H_

class Particle;

class Material_{
	
	protected:
	Particle *particle;
	public:
	Material(){}
	virtual inline double CalcYieldStress();
	virtual inline double CalcYieldStress(const double &strain, const double &strain_rate, const double &temp){};	
};

class JohnsonCook:
public Material_{
	double T_t,T_m;	//transition and melting temps
	double A, B, C;
	double n, m;
	double eps_0;
	
	public:
	JohnsonCook(){}
	//You provide the values of A, B, n, m, 
	//θmelt, and  θ_transition
	//as part of the metal plasticity material definition.
	JohnsonCook(const double &a, const double &b, const double &c, const double &eps_0):
	A(a),B(b),C(c){}
	inline double CalcYieldStress(){}	
	inline double CalcYieldStress(const double &strain, const double &strain_rate, const double &temp);	
};
#include "Material.cpp"

#endif