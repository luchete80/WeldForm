
#include "Domain.h"

using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
       SPH::Domain	dom;

        dom.Dimension	= 3;
        dom.Nproc	= 4;
    	dom.Kernel_Set(Qubic_Spline);
    	dom.Scheme	= 1;

        double dx,h,rho,K,G,Cs,Fy;
    	double H,L,n;

    	H	= 1.;
    	n	= 20.;

    	rho	= 1000.0;
    	dx	= H / n;
    	h	= dx*1.2; //Very important
        Cs	= sqrt(K/rho);

        double timestep;

        cout<<"t  = "<<timestep<<endl;
        cout<<"Cs = "<<Cs<<endl;
        cout<<"K  = "<<K<<endl;
        cout<<"G  = "<<G<<endl;
        cout<<"Fy = "<<Fy<<endl;
    	//dom.GeneralAfter = & UserAcc;
      
      ///// IMPORTANT
      dom.DomMax(0) = H;
      dom.DomMin(0) = -H;
      dom.ReadXDMF			("fsw.hdf5");
      
		std::cout << "Particle Number: "<< dom.Particles.size() << endl;
     	double x;
			dom.gradKernelCorr = false;
    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		x = dom.Particles[a]->x(0);
			dom.Particles[a]->k_T			=	3000.;
			dom.Particles[a]->cp_T			=	1.;
			dom.Particles[a]->h_conv		= 100.0; //W/m2-K
			dom.Particles[a]->T_inf 		= 500.;
			dom.Particles[a]->T				= 20.0;			
    		if ( x < dx ) {
    			dom.Particles[a]->ID 			= 2;
    			dom.Particles[a]->Thermal_BC 	= TH_BC_CONVECTION;
				// cout << "Particle " << a << "is convection BC" <<endl;
			}
    	}

        timestep = (0.3*h*h*rho*dom.Particles[0]->cp_T/dom.Particles[0]->k_T);	
				timestep = 0.001;	
		cout << "Time Step: "<<timestep<<endl;
		//timestep=1.e-6;
		//0.3 rho cp h^2/k
    
    
    dom.WriteXDMF    ("test.hdf5");
	
	

        return 0;
}

MECHSYS_CATCH
