
#include "Domain.h"

using std::cout;
using std::endl;

#include <vector>
//Axisymmetric SPH simulation of elasto-plastic contact
//in the low velocity impact

int main(int argc, char **argv) try
{
  SPH::Domain	dom;
  
  int Dimension = 2;
  
  dom.Dimension	= 2; // In fact is one
  dom.Nproc	= 1;
  dom.Kernel_Set(Qubic_Spline);
  dom.Scheme	= 1;	//Mod Verlet

	double dx,h,rho;
	double L;

	L	= 1.0;		
  rho	= 1.0;
  dx = 0.1;
  h	= dx*1.2; //Very important

  dom.DomMax(0) = L;
  dom.DomMin(0) = 0.;
  dom.DomMax(1) = L;
  dom.DomMin(1) = 0.;
  
  Vec3_t p0 = Vec3_t ( 2.0, 0., 0.);
										
	//dom.AddCylinderLength(1, Vec3_t(0.,0.,-L/10.), R, L + 2.*L/10.,  dx/2., rho, h, false); 

// inline void Domain::AddBoxLength(int tag, Vec3_t const & V, double Lx, double Ly, double Lz, 
									// double r, double Density, double h, int type, int rotation, bool random, bool Fixed)
                  
  dom.AddBoxLength(1 ,Vec3_t(0.,0.,0.), 
                      L + dx/10.0 , L + dx/10.0 ,  0 , 
                      dx/2.0 ,rho, h, 1 , 0 , false, false );
      
  //cout << "Particle count: "<<dom.Particles.Size()<<endl;
  dom.WriteXDMF("maz");

  for (int i = 0; i<dom.Particles.Size();i++) 
    cout << "i" << i<< ", x: "<<dom.Particles[i]->x(0)<<endl;
    
    
  cout << "Cell Init"<<endl;
	dom.CellInitiate();
  cout << "SMPairs size: "<<dom.SMPairs[0].Size()<<endl;
  cout << "Generating list"<<endl;
	dom.ListGenerate();
  cout << "Searching for Nbs"<<endl;
	dom.MainNeighbourSearch();
  dom.SaveNeighbourData();
  cout << "Done"<<endl;
  
  cout << "Done"<<endl;
  
  std::vector<double>  fx(dom.Particles.Size());
  std::vector<double> dfx(dom.Particles.Size());

  std::vector<double> dfx_c(dom.Particles.Size());
	
  std::vector<double>  gx(dom.Particles.Size());
  std::vector<double>  in(dom.Particles.Size());

 cout << "Calculating Corrected Kernel Gradient..."<<endl;	
	dom.CalcGradCorrMatrix();

  cout << "Calculating Kernel..."<<endl;
	for (int k=0; k<dom.Nproc;k++) {
    double mi,mj;
    double di,dj;
    SPH::Particle *P1,*P2;
    Vec3_t xij;
		for (size_t a=0; a<dom.SMPairs[k].Size();a++) {//Same Material Pairs, Similar to Domain::LastComputeAcceleration ()
			//cout << "a: " << a << "p1: " << SMPairs[k][a].first << ", p2: "<< SMPairs[k][a].second<<endl;
      int i = dom.SMPairs[k][a].first;
      int j = dom.SMPairs[k][a].second;
			P1	= dom.Particles[i];
			P2	= dom.Particles[j];
			xij	= P1->x - P2->x;
						
			mi = P1->Mass;
			mj = P2->Mass;	
      
      di = P1->Density;
			dj = P2->Density;	
      double K	= SPH::Kernel(Dimension, 0, norm(xij)/h, h);
      double GK = SPH::GradKernel(Dimension, 0, norm(xij)/h, h);
      
      cout << "r, K, GK: "<<norm(xij)/h<<", "<<K<<", "<<GK<<endl;
      cout <<"Vi"<<mj/dj<<endl;
      //mj/dj = dx * dx
      in[i] += mj/dj /** (1 + P2->x(0))*(1.+P2->x(1)) */ * K;
      in[j] += mi/di /** (1 + P1->x(0))*(1.+P1->x(1))*/ * K;
      
      fx[i] += mj/dj * (1 + P2->x(0))*(1.+P2->x(1)) * K;
      fx[j] += mi/di * (1 + P1->x(0))*(1.+P1->x(1)) * K;
      
      // gx[i] += /*mj/dj*/ /*dx * */P2->x(0)* K;
      // gx[j] += /*mi/di*/ /*dx * */P1->x(0)* K;
      
      dfx[i] += mj/dj * (1 + P2->x(0))*(1.+P2->x(1)) * GK * xij(0);
      dfx[j] -= mi/di * (1 + P1->x(0))*(1.+P1->x(1)) * GK * xij(0);
			
			Vec3_t GK_ci, GK_cj; 
			Mult(xij, GK * dom.Particles[i]->gradCorrM,GK_ci);
			Mult(xij, GK * dom.Particles[j]->gradCorrM,GK_cj);	

      dfx_c[i] += mj/dj * (1 + P2->x(0))*(1.+P2->x(1)) * GK_cj(0);
      dfx_c[j] -= mi/di * (1 + P1->x(0))*(1.+P1->x(1)) * GK_ci(0);
      
		} //Nproc //Pairs  
  }

	for (int i = 0; i<dom.Particles.Size();i++) {
		double x = dom.Particles[i]->x(0);
		double y = dom.Particles[i]->x(1);
		double K	= SPH::Kernel(Dimension, 0, 0, h);
		
		fx[i] += /*mj/dj */dx * dx * (1.+x)*(1.+y) * K;

	}
	
  cout << "Done."<<endl;
  //dom.m_kernel = SPH::iKernel(dom.Dimension,h);	

   cout << "i, x,y, anal, num, int, nb, "<< endl;  
  for (int i = 0; i<dom.Particles.Size();i++) {
    double x = dom.Particles[i]->x(0);
    double y = dom.Particles[i]->x(1);
    cout << i<<", "<<x<<", "<<y<<", "<<(1.+x)*(1.+y)<<", "<<fx[i]<<", "<<in[i]<<dom.Particles[i]->Nb<<endl;
  }

 cout << endl<< "i, x,y, grad anal, grad num, grad corr num, nb, "<< endl;  
  for (int i = 0; i<dom.Particles.Size();i++) {
    double x = dom.Particles[i]->x(0);
    double y = dom.Particles[i]->x(1);
    cout << i<<", "<<x<<", "<<y<<", "<<(1.+y)<<", "<<dfx[i]<<", "<<dfx_c[i]<<", "<<dom.Particles[i]->Nb<<endl;
  }
	
  // cout << endl<< "Derivatives"<<endl;
  // for (int i = 0; i<dom.Particles.Size();i++) {
    // double x = dom.Particles[i]->x(0);
    // cout << "Analytical" << dom.Particles[i]->x(0)<<", "<<x*x<<endl;
    // cout << dom.Particles[i]->x(0)<<", "<<dfx[i]<<endl;
  // }
  
  return 0;
}
MECHSYS_CATCH
