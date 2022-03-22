
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
  dx = 0.05;
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
  
  std::vector<Vec3_t>  vx(dom.Particles.Size());
  std::vector<Vec3_t> dfx(dom.Particles.Size());

  std::vector<Vec3_t> dfx_c(dom.Particles.Size());

  std::vector< Mat3_t > grad_vx(dom.Particles.Size());
  std::vector< Mat3_t > grad_vx_c(dom.Particles.Size());
	
  std::vector<double>  gx(dom.Particles.Size());
  
  std::vector<Vec3_t>  vij(dom.Particles.Size());

///////////////////////////////////////
///// DEFINING VELOCITY FIELD:

for (int i = 0; i<dom.Particles.Size();i++) {
		double x = dom.Particles[i]->x(0);
		double y = dom.Particles[i]->x(1);

		vx[i](0) = (y-0.5)*(y-0.5) ;
		vx[i](1) = x;
		vx[i](2) = 0;
	}

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

      //grad va = Sum b (vb-va) X gradWb(xa)
      //Bonet et. al 1999 eqn (43)
      Vec3_t vij = vx[i] - vx[j];
      Mat3_t t;
			Dyad (vij, Vec3_t(GK*xij),t);
      
      grad_vx[i] = grad_vx[i] + mj/dj * t;
      grad_vx[j] = grad_vx[j] - mi/di * t;
			
			Vec3_t GK_ci, GK_cj; 

			Mult(xij, GK * dom.Particles[i]->gradCorrM,GK_ci);
			Mult(xij, GK * dom.Particles[j]->gradCorrM,GK_cj);	
			Mat3_t t_ci,t_cj;
			Dyad (vij, GK_ci,t_ci);
			Dyad (vij, GK_cj,t_cj);
			//cout << "x1,x2,

			//This is like gradf_i = COrr x gradf
      grad_vx_c[i] = grad_vx_c[i] + mj/dj * t_ci; //If cj fails if ci is like applying 
      grad_vx_c[j] = grad_vx_c[j] - mi/di * t_cj;
			

			// Dyad (vxij, Vec3_t(GK*xij),t);
			// test [i] = test[i] - mj/dj * t;
			// test [j] = test[i] + mi/di * t;
      
		} //Nproc //Pairs  
  }

	for (int i = 0; i<dom.Particles.Size();i++) {
		double x = dom.Particles[i]->x(0);
		double y = dom.Particles[i]->x(1);
		double K	= SPH::Kernel(Dimension, 0, 0, h);
		
		//fx[i] += /*mj/dj */dx * dx * (1.+x)*(1.+y) * K;

	}
	
  cout << "Done."<<endl;
  //dom.m_kernel = SPH::iKernel(dom.Dimension,h);	

		// vx[i](0) = (y-0.5)*(y-0.5) + x;
		// vx[i](1) = x;
		
   cout << "i, x,y, dvxdy anal, num 01 num 10, nb, "<< endl;  
  for (int i = 0; i<dom.Particles.Size();i++) {
    double x = dom.Particles[i]->x(0);
    double y = dom.Particles[i]->x(1);
		
		// double GK = (1.+y)* SPH::GradKernel(Dimension, 0, 0., h);
		// Vec3_t GK_c; 
		// Mult(dom.Particles[i]->gradCorrM,dfx[i],GK_c);
		
    cout << i<<", "<<x<<", "<<y<<", "<< 2*(y-0.5)<<", "<<grad_vx[i](0,1)<<", "<<grad_vx[i](1,0)<<", "<<grad_vx_c[i](0,1)<<endl;
  }

 // cout << endl<< "i, x,y, grad anal, grad num, grad corr num, nb, "<< endl;  
  // for (int i = 0; i<dom.Particles.Size();i++) {
    // double x = dom.Particles[i]->x(0);
    // double y = dom.Particles[i]->x(1);
		
		// double GK = (1.+y)* SPH::GradKernel(Dimension, 0, 0., h);
		// Vec3_t GK_c; 
		// Mult(dom.Particles[i]->gradCorrM,dfx[i],GK_c);
			
    // cout << i<<", "<<x<<", "<<y<<", "<<(1.+y)<<", "<<dfx[i](0)<<", "<<GK_c(0)<<", "<<dfx_c[i](0)<<", "<<dom.Particles[i]->Nb<<endl;
  // }
	

  
  return 0;
}
MECHSYS_CATCH
