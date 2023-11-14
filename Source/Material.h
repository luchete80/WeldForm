#ifndef _MATERIAL_H_
#define _MATERIAL_H_
#include <string>

class Particle;

class Material_;
class JohnsonCook;

using namespace std;

// TODO: SEE WHAT ABOUT SMAX
//SOULD BE PASSED TO INHERITED RANKINE OR NOT?
//Johson Cook Damage Ratio according to Islam 2017 and 
class DamageModel {
protected:
  std::string m_dam_crit_str;
public:
  double Gf; //Fracture Energ
	double sigma_max;
	double delta_max;  // being calculated for example delta_max = 2GF/(sigmamax) 

  
  DamageModel(){};
  DamageModel(const double &smax_, const double &Gf_)
	:sigma_max(smax_),Gf(Gf_){}  
  //std::string & getDamageCriterion() {return m_dam_crit_str;}  
	virtual std::string getDamageCriterion() {return string("");}  
  virtual double CalcFractureStrain(const double &pl_strain){}
  virtual double CalcFractureStrain(const double &str_rate, const double &sig_as,const double &T){}
	virtual ~DamageModel(){}
  
};

class RankineDamage:
public DamageModel {
	public:
  RankineDamage(const double &smax_, const double &Gf_)
	:DamageModel(smax_,Gf_){
    m_dam_crit_str = "Rankine";
  }
	std::string getDamageCriterion() {return string("Rankine");} 
	double CalcFractureStrain(const double &pl_strain){}
	virtual ~RankineDamage(){}
};

class JohnsonCookDamage:
public DamageModel {
  double m_D1,m_D2,m_D3,m_D4,m_D5; 
  double m_eps0;  //REFERENCE STRAIN RATE, SAME AS JOHNSON COOK COUNTERPART, REDUNDANT
  JohnsonCook *mat; //ONLY TO OBTAIN EPS_0, NOT USED
	public:
  JohnsonCookDamage(const double &D1, const double &D2, const double &D3, const double &D4, const double &D5, const double &m_eps0_)
	:m_D1(D1),m_D2(D2),m_D3(D3),m_D4(D4),m_D5(D5), m_eps0(m_eps0_){
    m_dam_crit_str = "JohnsonCook";
  }
	std::string getDamageCriterion() {return string("JohnsonCook");} 
	double CalcFractureStrain(const double &str_rate, const double &sig_as,const double &T) { //sig_as is stress triaxiality
    return ( (m_D1 + m_D2 *exp(m_D3*sig_as)) * pow (1.0 + (str_rate/m_eps0), m_D4) * (1.0 + m_D5 * T)); //Islam 2017 eq. 35, 
  }
	virtual ~JohnsonCookDamage(){}
};

class Elastic_{
	private:
	double E_m, nu_m;	//Poisson and young
	double K_m, G_m;
	
	public:
	Elastic_(){}
	Elastic_(const double &e, const double &nu):E_m(e),nu_m(nu){}
	const double& E()const{return E_m;}
	
};

class Material_{
	
	protected:
	Elastic_ elastic_m;

	double E_m, nu;	//TODO, move to elastic class
	public:
	DamageModel *damage;
	Material_(){}
	Material_(const Elastic_ el):elastic_m(el){}
	virtual inline double CalcTangentModulus(){return 0.0;};
	virtual inline double CalcTangentModulus(const double &strain, const double &strain_rate, const double &temp){return 0.0;};
	virtual inline double CalcTangentModulus(const double &strain){return 0.0;};
	virtual inline double CalcYieldStress(){return 0.0;}
	virtual inline double CalcYieldStress(const double &strain){return 0.0;};
	virtual inline double CalcYieldStress(const double &strain, const double &strain_rate, const double &temp){return 0.0;}
  virtual double &getRefStrainRate() {}//only for JC
	const Elastic_& Elastic()const{return elastic_m;}
  //~Material_();
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
	double eps_0; //ONLY FOR JC DAMAGE , CORRECT THIS
	
	public:
	JohnsonCook(){}
	//You provide the values of A, B, n, m, 
	//θmelt, and  θ_transition
	//as part of the metal plasticity material definition.
	JohnsonCook(const Elastic_ &el,const double &a, const double &b, const double &n_, 
              const double &c, const double &eps_0_,
              const double &m_, const double &T_m_, const double &T_t_):
	Material_(el),A(a),B(b),C(c),
  m(m_),n(n_),eps_0(eps_0_),T_m(T_m_),T_t(T_t_)
  {}
	inline double CalcYieldStress(){return 0.0;}	
	inline double CalcYieldStress(const double &plstrain){
     double Et =0.;

    if (plstrain > 0.)
      Et = n * B * pow(plstrain,n-1.);
    else 
      Et = Elastic().E()*0.1; //ARBITRARY! TODO: CHECK MATHEMATICALLY
    return Et;
  } //TODO: SEE IF INCLUDE	
	inline double CalcYieldStress(const double &strain, const double &strain_rate, const double &temp);	
	inline double CalcTangentModulus(const double &strain, const double &strain_rate, const double &temp);
  double &getRefStrainRate(){return eps_0;}//only for JC
  //~JohnsonCook(){}
};

class Hollomon:
public Material_{
	double K, m;
	double eps0;
  double eps1;  //if has a perfectly plastic algoritm
  double sy0;
	
	public:
	Hollomon(){}
	//You provide the values of A, B, n, m, 
	//θmelt, and  θ_transition
	//as part of the metal plasticity material definition.
	//ASSUMING AT FIRST COEFFICIENTS ARE GIVEN TO TOTAL STRAIN-STRESS
	Hollomon(const double eps0_, const double &k_, const double &m_):
	K(k_), m(m_){ eps0 = eps0_;}
	Hollomon(const Elastic_ &el, const double sy0_, const double &k_, const double &m_);
  
	inline double CalcTangentModulus(const double &strain);
	inline double CalcYieldStress(){return 0.0;}	
	inline double CalcYieldStress(const double &strain);	
};

// class Bilinear_Hollomon:
// public Material_{
	// double K, m;
	// double eps0;
  // double eps1;  //if has a perfectly plastic algoritm
	
	// public:
	// Hollomon(){}
	// //You provide the values of A, B, n, m, 
	// //θmelt, and  θ_transition
	// //as part of the metal plasticity material definition.
	// //ASSUMING AT FIRST COEFFICIENTS ARE GIVEN TO TOTAL STRAIN-STRESS
	// Hollomon(const double eps0_, const double &k_, const double &m_):
	// K(k_), m(m_){ eps0 = eps0_;}
	// Hollomon(const Elastic_ &el, const double eps0_, const double &k_, const double &m_):
	// Material_(el),K(k_), m(m_){ eps0 = eps0_;}
	// inline double CalcTangentModulus(const double &strain);
	// inline double CalcYieldStress(){}	
	// inline double CalcYieldStress(const double &strain);	
// };

#include "Material.cpp"

#endif