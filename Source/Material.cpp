#include "Particle.h"

///////////////////////////////////////////////////
//sy = [A + B(epl^n)] [1 + C ln(e_dot pl/e_dot 0) (1 - pow)]

inline double JohnsonCook::CalcYieldStress(const double &strain, const double &strain_rate, const double &temp)	{
	double T_h = (temp - T_t) / (T_m - T_t);
	double sr = strain_rate;
	if (strain_rate == 0.0)
		sr = 0.001;
	
	double sy = (A+B*pow(strain, n))*(1.0 + C * log (sr/ eps_0) ) * (1.0 - pow(T_h,m));
	
	return sy;
}	

#include<iostream>
using namespace std;
inline double JohnsonCook::CalcTangentModulus(const double &plstrain, const double &strain_rate, const double &temp)	{
	double sy, T_h;
  //cout << "n, B, C, eps_0, T_t, m"<<n<<", "<<B<<", "<<C<<"eps0, "<<eps_0<<", "<<", "<<T_t<<", "<<m<<endl;
	T_h = (temp - T_t) / (T_m - T_t);
	
  //double sy = (A+B*pow(strain, n))*(1.0 + C * log (strain_rate/ eps_0) ) * (1.0 - pow(T_h,m));
  double Et =0.;

  if (plstrain > 0.)
    Et = n * B * pow(plstrain,n-1.)*(1.0 + C*log(strain_rate/ eps_0)) * (1.0-pow (T_h,m));
  else 
    Et = Elastic().E()*0.1; //ARBITRARY! TODO: CHECK MATHEMATICALLY
  return Et;
}	

inline double Hollomon::CalcYieldStress(const double &strain)	{
	double sy = K*pow(strain + eps0, m);
	return sy;
}	

inline double Hollomon::CalcTangentModulus(const double &strain) {
	double Et = K*m*pow(strain + eps0, (m-1.0));
	//cout << "ET: "<<Et<<endl;
	return Et;
}