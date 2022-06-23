#include "Particle.h"

///////////////////////////////////////////////////
//sy = [A + B(epl^n)] [1 + C ln(e_dot pl/e_dot 0) (1 - pow)]

inline double JohnsonCook::CalcYieldStress(const double &strain, const double &strain_rate, const double &temp)	{
	double T_h = (temp - T_t) / (T_m - T_t);
	double sr = strain_rate;
	if (strain_rate == 0.0)
		sr = 1.e-5;
	
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
//Case with plastic plateau 
Hollomon::Hollomon(const Elastic_ &el, const double sy0_, const double &k_, const double &m_):
Material_(el),K(k_), m(m_) {
  eps0 = sy0_/el.E(); 
  sy0  = sy0_;
  eps1 = pow(sy0_/k_, 1./m);
  cout << "eps_0 "<<eps0 << ", eps_1 "<<eps1<<endl;
  if (eps0 > eps1){
    throw new Fatal("ERROR, Hollomon material bad definition, please correct Yield Stress, Elastic Modulus or Material hardening constants.");
  }
}  
  
inline double Hollomon::CalcYieldStress(const double &strain)	{
  double sy;
  if (strain + eps0 > eps1) sy = K*pow(strain + eps0, m); //plateau surpassed. If no plateau, eps1=eps0 so 
  else                      sy = sy0; 
	return sy;
}	

inline double Hollomon::CalcTangentModulus(const double &strain) {
	double Et;
  if (strain + eps0 > eps1) Et = K*m*pow(strain + eps0, (m-1.0));
  else                      Et = 0.;
	//cout << "ET: "<<Et<<endl;
	return Et;
}