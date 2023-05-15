/***********************************************************************************
* PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids        * 
*             and soils) using Smoothed Particle Hydrodynamics method              *   
* Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo-Torres              *
*                                                                                  *
* This file is part of PersianSPH                                                  *
*                                                                                  *
* This is free software; you can redistribute it and/or modify it under the        *
* terms of the GNU General Public License as published by the Free Software        *
* Foundation; either version 3 of the License, or (at your option) any later       *
* version.                                                                         *
*                                                                                  *
* This program is distributed in the hope that it will be useful, but WITHOUT ANY  *
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  *
* PARTICULAR PURPOSE. See the GNU General Public License for more details.         *
*                                                                                  *
* You should have received a copy of the GNU General Public License along with     *
* PersianSPH; if not, see <http://www.gnu.org/licenses/>                           *
************************************************************************************/
#define ARC_RADIUS		0.004
#define WELD_LENGTH		0.070	//This include 2 * radiuos
#define PROBE_LENGTH	0.12
#define THICKNESS		0.006
#define HEAT_EFF		0.65
#define WELD_EGY		40000.0 //ES LA TOTAL (SIN SIMETRIA)
#define T_INI				90.0
#define T_AMB				20.0
#define T_WELD			25.0

#define WALL_HEIGHT		0.006

#define	WELDING_SPEED	(WELD_LENGTH-2.0*ARC_RADIUS)/(T_WELD/60.) 	//m/min, Son 55mm en 25 minutos

#include <fstream>      // std::ofstream

#include "Domain.h"

double last_time_saved = 0.;
Vec3_t meas_pos[3]; //Temperature measurement points
int meas_part_idx[3];

int meas_mov_part_idx[3];

bool endheat;

double rpart;

int find_part_idx(SPH::Domain & dom, Vec3_t pos){
	int ret;
	double mindist = 1.;
	#pragma omp parallel for schedule (static) num_threads(dom.Nproc)
	for (size_t a=0; a<dom.Particles.Size(); a++) {
		Vec3_t dist =pos - dom.Particles[a]->x;
		if ( norm(dist)<mindist ) {
			mindist=norm(dist);
			ret = a;
		}
		//if ()
	}
	return ret;
}

std::ofstream ofs,ofsmov;

void UserAcc(SPH::Domain & dom)
{
	int meas_mov_part_idx;
	
	double xstart = -WELD_LENGTH/2. + ARC_RADIUS;
	double x,y,z;

	double total_heatflux = HEAT_EFF * ( WELD_EGY/2.) / T_WELD;	//100kW
	//cout << "Total heat, t_heat"<<total_heatflux<<", "<< t_heat << endl;
	int heatflux_partcount = 0;
		
	double xsource = xstart + WELDING_SPEED/60. * dom.getTime();
	
	if ( dom.getTime() < T_WELD ) {
		#pragma omp parallel for schedule (static) num_threads(dom.Nproc)
		for (size_t i=0; i<dom.Particles.Size(); i++){
			x = dom.Particles[i]->x(0);
			y = dom.Particles[i]->x(1);
			z = dom.Particles[i]->x(2);
			//cout << "z "<< dom.Particles[i]->x(2)<<endl;
			if (y <= (ARC_RADIUS +ARC_RADIUS/20 )&& y >= -ARC_RADIUS && z == THICKNESS/2. -rpart + WALL_HEIGHT){
				// cout << "z "<< dom.Particles[i]->x(2)<<endl;
				//cout << "x source "<<xsource<<", x: "<<dom.Particles[i]->x(0)  << endl;
				dom.Particles[i]->q_source = 0.;
				dom.Particles[i]->ID = 1;
				double r = sqrt ( (x-xsource)*(x-xsource)+y*y); //ysource=0
				//if ( x > xsource - ARC_RADIUS && x < xsource + ARC_RADIUS){
				if ( r <= (ARC_RADIUS +ARC_RADIUS/20) ){
					dom.Particles[i]->ID = 3;
					heatflux_partcount++;
				}
			}
		}
	} else { //Si no hago esto la fuente sigue avanzando 
		if (!endheat){
			#pragma omp parallel for schedule (static) num_threads(dom.Nproc)
			for (size_t i=0; i<dom.Particles.Size(); i++){
				if (dom.Particles[i]->ID == 3){
					dom.Particles[i]->q_source = 0.;
					dom.Particles[i]->ID = 1;
				}
			}
			endheat = true;
		}
	}
	double source = total_heatflux/heatflux_partcount ; //surface=1m2
	//cout << "Heat source particle count: "<<heatflux_partcount<<endl;
		
	#pragma omp parallel for schedule (static) num_threads(dom.Nproc)
	for (size_t a=0; a<dom.Particles.Size(); a++)
		if (dom.Particles[a]->ID == 3)
			dom.Particles[a]->q_source = source * dom.Particles[a]->Density / dom.Particles[a]->Mass;	
		


	
	if ( dom.getTime() > last_time_saved ){
		int p = meas_part_idx[0];
		ofs << dom.getTime() << ", " <<dom.Particles[p]->T<<endl;
		cout << "Particle temp: "<<dom.Particles[p]->T<<endl;
		Vec3_t pos; 
		if (!endheat)
			pos = Vec3_t(xsource,0.015,0.0);
		else
			pos = Vec3_t(WELD_LENGTH/2. - ARC_RADIUS,0.015,0.0);
		p = find_part_idx (dom,pos);
		ofsmov << dom.getTime() << ", "<< dom.Particles[p]->x(0)
								<< ", "<< dom.Particles[p]->x(1) 
								<< ", "<< dom.Particles[p]->x(2) 
								<< ", "<< dom.Particles[p]->T<<endl;	
		
		last_time_saved += 1.;
	}}


using std::cout;
using std::endl;

int main(int argc, char **argv) try
{
       SPH::Domain	dom;

        dom.Dimension	= 3;
        dom.Nproc	= 4;
    	dom.Kernel_Set(Qubic_Spline);
    	dom.Scheme	= 0;
//     	dom.XSPH	= 0.5; //Very important

        double dx,h,rho,K,G,Cs,Fy;
    	double W,L,T,n;//
		double th_factor = 0.5;

    	W	= 0.025;
    	L	= PROBE_LENGTH;
		T	= THICKNESS;
		n 	= 8.;
		
    	rho	= 7850.0;
    	dx	= T/n;

    	h	= dx*1.1; //Very important
        Cs	= sqrt(K/rho);

        double timestep;

        cout<<"t  = "<<timestep<<endl;
        cout<<"Cs = "<<Cs<<endl;
        cout<<"K  = "<<K<<endl;
        cout<<"G  = "<<G<<endl;
        cout<<"Fy = "<<Fy<<endl;
    	dom.GeneralAfter = & UserAcc;
        dom.DomMax(0) = L;
        dom.DomMin(0) = -L;
		
		rpart = dx/2.;
		

		if (WALL_HEIGHT == 0.0){
     	dom.AddBoxLength(1 ,Vec3_t ( -L/2.0 , 0.0 , -T/2.0 ), 
										L +dx, W +dx,  T +dx, 
										dx/2.0 , rho, h, 1 , 0 , false, false );
		
		} else {
			dom.Add3DCubicBoxParticles(1, Vec3_t ( -L/2.0 /*-dx/2.*/, 0.0 /*-dx/2.*/, -T/2.0 /*-dx/2.*/ ),
										L , W ,  T ,  
										dx/2., rho, h);

			int idx = find_part_idx(dom , Vec3_t(-WELD_LENGTH/2.0,0,-T/2.0));
			cout << "Init Partcile: "<< idx<< ", "<<dom.Particles[idx]->x(0)<<endl;
			Vec3_t pos_ini = Vec3_t (dom.Particles[idx]->x(0)-dx/2.,0.0,T/2.);
			
			dom.Add3DCubicBoxParticles(1, pos_ini,
										WELD_LENGTH +dx, ARC_RADIUS,  WALL_HEIGHT,  
										dx/2., rho, h);

			dom.Calculate3DMass(rho);									
		}
		std::cout << "Particle Number: "<< dom.Particles.size() << endl;
     	double x,y,z;
		int conv_partcount = 0;
    	for (size_t a=0; a<dom.Particles.Size(); a++)
    	{
    		x = dom.Particles[a]->x(0);
    		y = dom.Particles[a]->x(1);
    		z = dom.Particles[a]->x(2);
			dom.Particles[a]->k_T			=	50.;
			dom.Particles[a]->cp_T			=	490.;
			dom.Particles[a]->h_conv		= 	20.0; //W/m2-K
			dom.Particles[a]->T_inf 		= 	T_AMB;
			dom.Particles[a]->T				= 	T_INI;			
    		dom.Particles[a]->PresEq	= 0;
    		dom.Particles[a]->Cs		= Cs;
    		dom.Particles[a]->Shepard	= false;
    		dom.Particles[a]->Alpha		= 0.0;
    		dom.Particles[a]->Beta		= 0.0;
			
			//cout << "z "<< dom.Particles[a]->x(2)<<endl;
			double r=dx/2.;
    		if ( z == -T/2.+ r //Bottom face
				|| ( ((T/2.0 - r) - z )<1e-5 && y > ARC_RADIUS  )  //Top Face
				|| abs((W) - y) < dx //Side face
				|| x == (-L/2. + r ) 
				|| abs(L/2.-r-x)<dx 
				|| ( abs(y - ARC_RADIUS) < dx && z > T/2.0 && z < T/2 -r + WALL_HEIGHT)
				|| ( ((T/2.0 - r) -z )<1e-5 && abs(x)> WELD_LENGTH/2. ) 
				|| (  z > T/2.0  && ( abs(x) > (WELD_LENGTH/2. - r ) ) )  ){
    			dom.Particles[a]->ID 			= 2;
    			dom.Particles[a]->Thermal_BC 	= TH_BC_CONVECTION;
				// cout << "Particle " << a << "is convection BC" <<endl;
				conv_partcount++;
				}
				
				// if ( abs((W) - y) < dx){
					// cout << "y: "<<y<<endl;
				// }

    	}
		meas_pos[0](0)=0.0; meas_pos[0](1)=0.015;meas_pos[0](2)=0.0;
		
		for (int p=0;p<1;p++){
			double mindist=1.;
			for (size_t a=0; a<dom.Particles.Size(); a++) {
				
				
				Vec3_t dist =meas_pos[p] - dom.Particles[a]->x;
				if ( norm(dist)<mindist ) {
					mindist=norm(dist);
					meas_part_idx[p] = a;
				}
				//if ()
			}
		}
		cout << " Particle measured is "<<meas_part_idx[0]<<endl;
		int a = meas_part_idx[0];
		cout << "x,y,<:" <<dom.Particles[a]->x(0)<<", "<<dom.Particles[a]->x(1)<<", "<<dom.Particles[a]->x(2)<<endl;
		
		cout << "Convection particle count: "<<conv_partcount<<endl;
		
		
        timestep = th_factor * (0.3*h*h*rho*dom.Particles[0]->cp_T/dom.Particles[0]->k_T);	
		cout << "Time Step: "<<timestep<<endl;
		//timestep=1.e-6;
		//0.3 rho cp h^2/k
	
		
//    	dom.WriteXDMF("maz");
//    	dom.Solve(/*tf*/0.01,/*dt*/timestep,/*dtOut*/0.001,"test06",999);

		//dom.ThermalSolve(/*tf*/1.01,/*dt*/timestep,/*dtOut*/0.1,"test06",999);
		
		ofs.open   ("test.txt", std::ofstream::out | std::ofstream::app);
		ofsmov.open("test_mov.txt", std::ofstream::out | std::ofstream::app);
		
		endheat = false;
		
		dom.ThermalSolve(/*tf*/120.1,/*dt*/timestep,/*dtOut*/1.0,"test06",999);

		ofs.close();
		ofsmov.close();
        return 0;
}


			// dom.Particles[a]->k_T			=3000.;
			// dom.Particles[a]->cp_T			=1.;
			// dom.Particles[a]->h_conv		= 100.0; //W/m2-K
			// dom.Particles[a]->T_inf 		= 500.;
			// dom.Particles[a]->T				= 20.0;
    		// x = dom.Particles[a]->x(0);
    		// if (x=-H/2.0)
    			// dom.Particles[a]->Thermal_BC=TH_BC_CONVECTION;

MECHSYS_CATCH
