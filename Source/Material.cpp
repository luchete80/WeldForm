#include "Particle.h"

///////////////////////////////////////////////////
//sy = [A + B(epl^n)] [1 + C ln(e_dot pl/e_dot 0) (1 - pow)]

inline double JohnsonCook::CalcYieldStress(const double &strain, const double &strain_rate, const double &temp)	{
  // OLD////////////////////////
	// double T_h = (temp - T_t) / (T_m - T_t);
	// double sr = strain_rate;
	// if (strain_rate == 0.0)
		// sr = 1.e-5;
	
	// double sy = (A+B*pow(strain, n))*(1.0 + C * log (sr/ eps_0) ) * (1.0 - pow(T_h,m));  
  
  // NEW /////////////////////
	double T_h = (temp - T_t) / (T_m - T_t);
	double sr = strain_rate;
  double f = 1.0;
	if (strain_rate > eps_0)
		f = (1.0 + C * log(strain_rate/eps_0));
	
	double sy = (A+B*pow(strain, n)) * f * (1.0 - pow(T_h,m));
	
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
  double f = 1.0;
	if (strain_rate > eps_0)
		f = (1.0 + C * log(strain_rate/eps_0));
   if (plstrain > 0.)
    Et = n * B * pow(plstrain,n-1.) * f * (1.0-pow (T_h,m));
   else 
     //Et = Elastic().E()*0.1; //ARBITRARY! TODO: CHECK MATHEMATICALLY
		Et = Elastic().E();
	 
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



inline double JohnsonCookDamage::CalcFractureStrain(const double &str_rate, const double &sig_as,const double &T) { //sig_as is stress triaxiality
	//return ( (m_D1 + m_D2 *exp(m_D3*sig_as)) * pow (1.0 + (str_rate/m_eps0), m_D4) * (1.0 + m_D5 * T)); //Islam 2017 eq. 35, 
	double T_h = (T - mat->getT_t()) / (mat->getT_m() - mat->getT_t());
	//cout << "T_t, T_m "<< mat->getT_t() << ", "<<mat->getT_m()<<endl;
	double f_sr = 1.0;
	double f_triax = m_D1;
	if (str_rate > m_eps0) f_sr = 1.0 + log(str_rate/m_eps0)* m_D4;
	if (sig_as<1.5) f_triax = m_D1 + m_D2 *exp(m_D3*sig_as); //Pressure positive
	double ef = f_triax * f_sr * (1.0 + m_D5 * T_h);
	// cout << "sig_as" << sig_as<<endl;
	// if (ef <1.0 )ef = 1.0;
	//cout << "ef " << ef<<endl;
	if (ef>10.0) ef = 10.0;
	return ef; //ABAQUS
}

///////////////////////////////////////////////////
//sy = [A + B(epl^n)] [1 + C ln(e_dot pl/e_dot 0) (1 - pow)]
////////////// TEMPERATURE SHOULD BE IN CELSIUS

inline double GMT::CalcYieldStress(const double &strain, const double &strain_rate, const double &temp)	{
  // OLD////////////////////////
	// double T_h = (temp - T_t) / (T_m - T_t);
	// double sr = strain_rate;
	// if (strain_rate == 0.0)
		// sr = 1.e-5;
    
	// double sy = (A+B*pow(strain, n))*(1.0 + C * log (sr/ eps_0) ) * (1.0 - pow(T_h,m));  
  
  // NEW /////////////////////
	double T_h = (temp - T_t) / (T_m - T_t);
	double sr = strain_rate;
  double f = 1.0;
	// if (strain_rate > eps_0)
		// f = (1.0 + C * log(strain_rate/eps_0));
	
	double e,er,T, sy;
  e = strain; er = strain_rate; T = temp;
  if      (strain < e_min) e = e_min;
  else if (strain > e_max) e = e_max;

  if      (strain_rate < er_min) er = er_min;
  else if (strain_rate > er_max) er = er_max;

  if      (temp < er_min) T = T_min;
  else if (temp > er_max) T = T_max;
  
  sy = C1 * exp(C2*T)*pow(e,n1*T+n2) * exp((I1*T+I2)/e)*pow(er,m1*T+m2);
	
	return sy;
}	

inline double GMT::CalcTangentModulus(const double &plstrain, const double &strain_rate, const double &temp)	{
	double T_h;
  
	double e,er,T, sy;
  e = plstrain; er = strain_rate; T = temp;
  if      (plstrain < e_min) e = e_min;
  else if (plstrain > e_max) e = e_max;

  if      (strain_rate < er_min) er = er_min;
  else if (strain_rate > er_max) er = er_max;

  if      (temp < er_min) T = T_min;
  else if (temp > er_max) T = T_max;
  
  //double sy = (A+B*pow(strain, n))*(1.0 + C * log (strain_rate/ eps_0) ) * (1.0 - pow(T_h,m));
  double Et =0.;

	Et = C1*exp(C2*T)*pow(er,m1*T+m2)* //constant part
       pow(e,T*n1+n2-2.0)*(-I1*T-I2+e*(n1*T+n2))*exp((I1*T+I2)/e);

  // if (strain_rate > eps_0)
		// f = (1.0 + C * log(strain_rate/eps_0));
   // if (plstrain > 0.)
    // Et = n * B * pow(plstrain,n-1.) * f * (1.0-pow (T_h,m));
   // else 
     // //Et = Elastic().E()*0.1; //ARBITRARY! TODO: CHECK MATHEMATICALLY
		// Et = Elastic().E();
	 
  return Et;
}	