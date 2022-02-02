#ifndef _MATERIAL_H_
#define _MATERIAL_H_

class Particle;

class Elastic_{
	private:
	double E_m, nu_m;	//Poisson and young
	double K_m, G_m;
	
	public:
	const double& E()const{return E_m;}
	
};

class Material_{
	
	protected:
	Elastic_ elastic_m;
	public:
	Material(){}
	virtual inline double CalcYieldStress();
	virtual inline double CalcYieldStress(const double &strain, const double &strain_rate, const double &temp){}
	const Elastic_& Elastic()const{return elastic_m;}
};

class _Plastic{
	
	public:
	virtual inline double CalcYieldStress();	
	//virtual inline double CalcYieldStress();
};

//TODO: derive johnson cook as plastic material flow
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

class Hollomon:
public Material_{
	double n, m;
	double eps_0;
	
	public:
	Hollomon(){}
	//You provide the values of A, B, n, m, 
	//θmelt, and  θ_transition
	//as part of the metal plasticity material definition.
	Hollomon(const double &m_):
	m(m_){}
	inline double CalcYieldStress(){}	
	inline double CalcYieldStress(const double &strain);	
};

#include "Material.cpp"

#endif