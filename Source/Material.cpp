#include "Particle.h"

///////////////////////////////////////////////////
//sy = [A + B(epl^n)] [1 + C ln(e_dot pl/e_dot 0) (1 - pow)]

inline double JohnsonCook::CalcYieldStress(const double &strain, const double &strain_rate, const double &temp)	{
	double T_h = (temp - T_t) / (T_m - T_t);
	double sy = (A+B*pow(strain, n))*(1.0 + C * log (strain_rate/ eps_0) ) * (1.0 - pow(T_h,m));
	return sy;
}	

inline double JohnsonCook::CalcTangentModulus(const double &strain, const double &strain_rate, const double &temp)	{
	double sy;
	//double T_h = (temp - T_t) / (T_m - T_t);
	//double sy = (A+B*pow(strain, n))*(1.0 + C * log (strain_rate/ eps_0) ) * (1.0 - pow(T_h,m));
	return sy;
}	

inline double Hollomon::CalcYieldStress(const double &strain)	{
	double sy = K*pow(strain, m);
	return sy;
}	

inline double Hollomon::CalcTangentModulus(const double &strain) {
	double Et = K*m*pow(strain + eps0, (m-1.0));
	//cout << "ET: "<<Et<<endl;
	return Et;
}