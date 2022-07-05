[1mdiff --git a/CMakeLists.txt b/CMakeLists.txt[m
[1mindex 7332705..ae1a175 100644[m
[1m--- a/CMakeLists.txt[m
[1m+++ b/CMakeLists.txt[m
[36m@@ -94,14 +94,19 @@[m [mSET(EXES[m
   #Contact_Compression_2surf[m
   #Contact_Compression_vel[m
   #Contact_Compression_test[m
[31m-  #Contact_Compression_half_ghost[m
[32m+[m[32m  Contact_Compression_half_ghost[m
   #Contact_Compression_half_ghost_debug[m
   #Contact_Compression_qt_ghost[m
   #Contact_Compression_TRIP[m
   #Contact_Compression_TRIP_half_ghost[m
   #Contact_Compression_TRIP_2surf[m
[32m+[m[32m  #Contact_Compression_1010_hlaf_ghost[m
[32m+[m[32m  #Compression_1010[m
[32m+[m[32m  #Contact_Compression_1010[m
[32m+[m[32m  #Compression_1010_half_ghost[m
 	#NastranReader_test[m
 	#CompressionForces[m
[32m+[m[32m  #Compression[m
 	#CompressionVel[m
   #CompressionVel_Implicit[m
 	#CompressionForces_Ghost[m
[36m@@ -126,7 +131,7 @@[m [mSET(EXES[m
   #Kernel_Test[m
 	#Kernel_Test_2D[m
 	#Kernel_Test_2D_vel[m
[31m-	Traction[m
[32m+[m	[32m#Traction[m
 	#traction_app[m
 	#Traction_Hollomon[m
 	#Traction_JC[m
[36m@@ -165,7 +170,9 @@[m [mFOREACH(var ${EXES})[m
     SET_TARGET_PROPERTIES (${var} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")[m
 ENDFOREACH(var)[m
 [m
[31m-	[m
[32m+[m[32mADD_EXECUTABLE        (WeldForm "${PROJECT_SOURCE_DIR}/Source/WeldForm.cpp" )[m
[32m+[m[32mTARGET_LINK_LIBRARIES (WeldForm ${LIBS} blitz)[m
[32m+[m[32mSET_TARGET_PROPERTIES (WeldForm PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")[m[41m	[m
     # ADD_EXECUTABLE        (${EXES} "Plate-Yield.cpp")[m
     # TARGET_LINK_LIBRARIES (${EXES} ${LIBS} blitz)[m
     # SET_TARGET_PROPERTIES (${EXES} PROPERTIES COMPILE_FLAGS "${FLAGS}" LINK_FLAGS "${LFLAGS}")[m
\ No newline at end of file[m
[1mdiff --git a/Compression.cpp b/Compression.cpp[m
[1mnew file mode 100644[m
[1mindex 0000000..35afcf0[m
[1m--- /dev/null[m
[1m+++ b/Compression.cpp[m
[36m@@ -0,0 +1,174 @@[m
[32m+[m
[32m+[m[32m/***********************************************************************************[m
[32m+[m[32m* PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids        *[m[41m [m
[32m+[m[32m*             and soils) using Smoothed Particle Hydrodynamics method              *[m[41m   [m
[32m+[m[32m* Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo-Torres              *[m
[32m+[m[32m*                                                                                  *[m
[32m+[m[32m* This file is part of PersianSPH                                                  *[m
[32m+[m[32m*                                                                                  *[m
[32m+[m[32m* This is free software; you can redistribute it and/or modify it under the        *[m
[32m+[m[32m* terms of the GNU General Public License as published by the Free Software        *[m
[32m+[m[32m* Foundation; either version 3 of the License, or (at your option) any later       *[m
[32m+[m[32m* version.                                                                         *[m
[32m+[m[32m*                                                                                  *[m
[32m+[m[32m* This program is distributed in the hope that it will be useful, but WITHOUT ANY  *[m
[32m+[m[32m* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  *[m
[32m+[m[32m* PARTICULAR PURPOSE. See the GNU General Public License for more details.         *[m
[32m+[m[32m*                                                                                  *[m
[32m+[m[32m* You should have received a copy of the GNU General Public License along with     *[m
[32m+[m[32m* PersianSPH; if not, see <http://www.gnu.org/licenses/>                           *[m
[32m+[m[32m************************************************************************************/[m
[32m+[m
[32m+[m[32m#include "Domain.h"[m
[32m+[m[32m#include "InteractionAlt.cpp"[m
[32m+[m
[32m+[m[32m#define TAU		0.005[m
[32m+[m[32m#define VMAX	10.0[m
[32m+[m
[32m+[m[32mint forcepart_count;[m
[32m+[m
[32m+[m
[32m+[m[32mvoid UserAcc(SPH::Domain & domi)[m
[32m+[m[32m{[m
[32m+[m	[32mdouble vcompress;[m
[32m+[m
[32m+[m
[32m+[m	[32mdouble acc = 7.906e4/(forcepart_count*domi.Particles[0]->Mass);[m
[32m+[m[41m	[m
[32m+[m	[32mif (domi.getTime() < TAU )[m[41m [m
[32m+[m		[32mvcompress = VMAX/TAU * domi.getTime();[m
[32m+[m	[32melse[m
[32m+[m		[32mvcompress = VMAX;[m
[32m+[m	[32m//cout << "time: "<< domi.getTime() << "V compress "<< vcompress <<endl;[m
[32m+[m	[32m#pragma omp parallel for schedule (static) num_threads(domi.Nproc)[m
[32m+[m
[32m+[m	[32m#ifdef __GNUC__[m
[32m+[m	[32mfor (size_t i=0; i<domi.Particles.Size(); i++)[m
[32m+[m	[32m#else[m
[32m+[m	[32mfor (int i=0; i<domi.Particles.Size(); i++)[m
[32m+[m	[32m#endif[m
[32m+[m[41m	[m
[32m+[m	[32m{[m
[32m+[m		[32mif (domi.Particles[i]->ID == 3)[m
[32m+[m		[32m{[m
[32m+[m			[32m//domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);[m
[32m+[m			[32m//domi.Particles[i]->a(2)		= -VMAX/TAU;[m
[32m+[m			[32m// domi.Particles[i]->v		= Vec3_t(0.0,0.0,-vcompress);[m
[32m+[m			[32m// domi.Particles[i]->va		= Vec3_t(0.0,0.0,-vcompress);[m
[32m+[m			[32m// domi.Particles[i]->vb		= Vec3_t(0.0,0.0,-vcompress);[m
[32m+[m			[32mdomi.Particles[i]->v(2)			= -vcompress;[m
[32m+[m			[32mdomi.Particles[i]->va(2)		= -vcompress;[m
[32m+[m			[32mdomi.Particles[i]->vb(2)		= -vcompress;[m
[32m+[m[32m//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);[m
[32m+[m		[32m}[m
[32m+[m		[32mif (domi.Particles[i]->ID == 2)[m
[32m+[m		[32m{[m
[32m+[m			[32mdomi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);[m
[32m+[m			[32mdomi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);[m
[32m+[m			[32mdomi.Particles[i]->vb		= Vec3_t(0.0,0.0,0.0);[m
[32m+[m			[32m//domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);[m
[32m+[m		[32m}[m
[32m+[m	[32m}[m
[32m+[m[32m}[m
[32m+[m
[32m+[m
[32m+[m[32musing std::cout;[m
[32m+[m[32musing std::endl;[m
[32m+[m
[32m+[m[32mint main(int argc, char **argv) try[m
[32m+[m[32m{[m
[32m+[m[32m  SPH::Domain	dom;[m
[32m+[m
[32m+[m[32m  dom.Dimension	= 3;[m
[32m+[m[32m  dom.Nproc	= 4;[m
[32m+[m[32m  dom.Kernel_Set(Qubic_Spline);[m
[32m+[m[32m  //dom.Kernel_Set(Hyperbolic_Spline);[m
[32m+[m[32m  dom.Scheme	= 1;	//Mod Verlet[m
[32m+[m[32m  //dom.XSPH	= 0.1; //Very important[m
[32m+[m
[32m+[m[32m  double dx,h,rho,K,G,Cs,Fy;[m
[32m+[m[32m  double R,L,n;[m
[32m+[m
[32m+[m[32m  R	= 0.15;[m
[32m+[m[32m  L	= 0.56;[m
[32m+[m[32m  n	= 30.0;		//in length, radius is same distance[m
[32m+[m
[32m+[m[32m  rho	= 2700.0;[m
[32m+[m[32m  K	= 6.7549e10;[m
[32m+[m[32m  G	= 2.5902e10;[m
[32m+[m[32m  Fy	= 300.e6;[m
[32m+[m[32m  //dx	= L / (n-1);[m
[32m+[m[32m  //dx = L/(n-1);[m
[32m+[m[32m  dx = 0.015;[m
[32m+[m[32m  h	= dx*1.2; //Very important[m
[32m+[m[32m    Cs	= sqrt(K/rho);[m
[32m+[m
[32m+[m[32m    double timestep;[m
[32m+[m[32m    timestep = (0.2*h/(Cs));[m
[32m+[m
[32m+[m[32m  //timestep = 2.5e-6;[m
[32m+[m
[32m+[m[32m  cout<<"t  = "<<timestep<<endl;[m
[32m+[m[32m  cout<<"Cs = "<<Cs<<endl;[m
[32m+[m[32m  cout<<"K  = "<<K<<endl;[m
[32m+[m[32m  cout<<"G  = "<<G<<endl;[m
[32m+[m[32m  cout<<"Fy = "<<Fy<<endl;[m
[32m+[m[32m  dom.GeneralAfter = & UserAcc;[m
[32m+[m[32m  dom.DomMax(0) = L;[m
[32m+[m[32m  dom.DomMin(0) = -L;[m
[32m+[m
[32m+[m
[32m+[m		[32m// inline void Domain::AddCylinderLength(int tag, Vec3_t const & V, double Rxy, double Lz,[m[41m [m
[32m+[m									[32m// double r, double Density, double h, bool Fixed) {[m
[32m+[m[41m										[m
[32m+[m		[32mdom.AddCylinderLength(1, Vec3_t(0.,0.,-L/10.), R, L + 2.*L/10.,  dx/2., rho, h, false);[m[41m [m
[32m+[m[41m		[m
[32m+[m		[32mcout << "Particle count: "<<dom.Particles.Size()<<endl;[m
[32m+[m[41m		[m
[32m+[m		[32mforcepart_count = 0;[m
[32m+[m		[32m//dom.gradKernelCorr = true;[m
[32m+[m		[32mdom.ts_nb_inc = 5;[m[41m		[m
[32m+[m[41m		[m
[32m+[m[41m    [m	[32mfor (size_t a=0; a<dom.Particles.Size(); a++)[m
[32m+[m[41m    [m	[32m{[m
[32m+[m[41m    [m		[32mdom.Particles[a]->G		= G;[m
[32m+[m[41m    [m		[32mdom.Particles[a]->PresEq	= 0;[m
[32m+[m[41m    [m		[32mdom.Particles[a]->Cs		= Cs;[m
[32m+[m[41m    [m		[32mdom.Particles[a]->Shepard	= false;[m
[32m+[m[41m    [m		[32mdom.Particles[a]->Material	= 2;[m
[32m+[m[41m    [m		[32mdom.Particles[a]->Fail		= 1;[m
[32m+[m[41m    [m		[32mdom.Particles[a]->Sigmay	= Fy;[m
[32m+[m[41m    [m		[32mdom.Particles[a]->Alpha		= 1.0;[m
[32m+[m[41m    [m		[32m//dom.Particles[a]->Beta		= 1.0;[m
[32m+[m[41m    [m		[32mdom.Particles[a]->TI		= 0.3;[m
[32m+[m[41m    [m		[32mdom.Particles[a]->TIInitDist	= dx;[m
[32m+[m[41m    [m		[32mdouble z = dom.Particles[a]->x(2);[m
[32m+[m[41m    [m		[32mif ( z < 0 ){[m
[32m+[m[41m    [m			[32mdom.Particles[a]->ID=2;[m
[32m+[m	[41m    [m			[32m// dom.Particles[a]->IsFree=false;[m
[32m+[m[41m    [m			[32m// dom.Particles[a]->NoSlip=true;[m[41m			[m
[32m+[m[41m				[m
[32m+[m				[32m}[m
[32m+[m[41m    [m		[32mif ( z > (L - dx) ) {//Changed to only last row[m
[32m+[m[41m    [m			[32mdom.Particles[a]->ID=3;[m
[32m+[m					[32m//dom.Particles[a]->XSPH		= 0.1;[m
[32m+[m					[32mforcepart_count++;[m
[32m+[m				[32m}[m
[32m+[m[41m    [m	[32m}[m
[32m+[m			[32mcout << "Contact Force Particles: "<<forcepart_count<<endl;[m
[32m+[m		[32mdom.WriteXDMF("maz");[m
[32m+[m		[32mdom.m_kernel = SPH::iKernel(dom.Dimension,h);[m[41m	[m
[32m+[m		[32mdom.BC.InOutFlow = 0;[m
[32m+[m
[32m+[m[32m    //dom.Solve_orig_Ext(/*tf*/0.00205,/*dt*/timestep,/*dtOut*/0.001,"test06",999);[m
[32m+[m		[32m//dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/0.0001,"test06",999);[m
[32m+[m[41m    [m
[32m+[m[32m    timestep = (0.2*h/(Cs+VMAX));[m
[32m+[m[41m    [m
[32m+[m[32m    dom.auto_ts = false;[m
[32m+[m[32m    //dom.Solve(/*tf*/0.0105,/*dt*/timestep,/*dtOut*/0.0001,"test06",999);[m
[32m+[m[32m    dom.SolveDiffUpdateKickDrift(/*tf*/0.105,/*dt*/timestep,/*dtOut*/1.e-4,"test06",10000);[m[41m	[m
[32m+[m[41m  [m
[32m+[m		[32mreturn 0;[m
[32m+[m[32m}[m
[32m+[m[32mMECHSYS_CATCH[m
[1mdiff --git a/CompressionVel.cpp b/CompressionVel.cpp[m
[1mnew file mode 100644[m
[1mindex 0000000..0e073da[m
[1m--- /dev/null[m
[1m+++ b/CompressionVel.cpp[m
[36m@@ -0,0 +1,172 @@[m
[32m+[m
[32m+[m[32m/***********************************************************************************[m
[32m+[m[32m* PersianSPH - A C++ library to simulate Mechanical Systems (solids, fluids        *[m[41m [m
[32m+[m[32m*             and soils) using Smoothed Particle Hydrodynamics method              *[m[41m   [m
[32m+[m[32m* Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo-Torres              *[m
[32m+[m[32m*                                                                                  *[m
[32m+[m[32m* This file is part of PersianSPH                                                  *[m
[32m+[m[32m*                                                                                  *[m
[32m+[m[32m* This is free software; you can redistribute it and/or modify it under the        *[m
[32m+[m[32m* terms of the GNU General Public License as published by the Free Software        *[m
[32m+[m[32m* Foundation; either version 3 of the License, or (at your option) any later       *[m
[32m+[m[32m* version.                                                                         *[m
[32m+[m[32m*                                                                                  *[m
[32m+[m[32m* This program is distributed in the hope that it will be useful, but WITHOUT ANY  *[m
[32m+[m[32m* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  *[m
[32m+[m[32m* PARTICULAR PURPOSE. See the GNU General Public License for more details.         *[m
[32m+[m[32m*                                                                                  *[m
[32m+[m[32m* You should have received a copy of the GNU General Public License along with     *[m
[32m+[m[32m* PersianSPH; if not, see <http://www.gnu.org/licenses/>                           *[m
[32m+[m[32m************************************************************************************/[m
[32m+[m
[32m+[m[32m#include "Domain.h"[m
[32m+[m[32m#include "InteractionAlt.cpp"[m
[32m+[m
[32m+[m[32m#define TAU		0.005[m
[32m+[m[32m#define VMAX	10.0[m
[32m+[m
[32m+[m[32mint forcepart_count;[m
[32m+[m
[32m+[m
[32m+[m[32mvoid UserAcc(SPH::Domain & domi)[m
[32m+[m[32m{[m
[32m+[m	[32mdouble vcompress;[m
[32m+[m
[32m+[m
[32m+[m	[32mdouble acc = 7.906e4/(forcepart_count*domi.Particles[0]->Mass);[m
[32m+[m[41m	[m
[32m+[m	[32mif (domi.getTime() < TAU )[m[41m [m
[32m+[m		[32mvcompress = VMAX/TAU * domi.getTime();[m
[32m+[m	[32melse[m
[32m+[m		[32mvcompress = VMAX;[m
[32m+[m	[32m//cout << "time: "<< domi.getTime() << "V compress "<< vcompress <<endl;[m
[32m+[m	[32m#pragma omp parallel for schedule (static) num_threads(domi.Nproc)[m
[32m+[m
[32m+[m	[32m#ifdef __GNUC__[m
[32m+[m	[32mfor (size_t i=0; i<domi.Particles.Size(); i++)[m
[32m+[m	[32m#else[m
[32m+[m	[32mfor (int i=0; i<domi.Particles.Size(); i++)[m
[32m+[m	[32m#endif[m
[32m+[m[41m	[m
[32m+[m	[32m{[m
[32m+[m		[32mif (domi.Particles[i]->ID == 3)[m
[32m+[m		[32m{[m
[32m+[m			[32m//domi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);[m
[32m+[m			[32m//domi.Particles[i]->a(2)		= -VMAX/TAU;[m
[32m+[m			[32m// domi.Particles[i]->v		= Vec3_t(0.0,0.0,-vcompress);[m
[32m+[m			[32m// domi.Particles[i]->va		= Vec3_t(0.0,0.0,-vcompress);[m
[32m+[m			[32m// domi.Particles[i]->vb		= Vec3_t(0.0,0.0,-vcompress);[m
[32m+[m			[32mdomi.Particles[i]->v(2)			= -vcompress;[m
[32m+[m			[32mdomi.Particles[i]->va(2)		= -vcompress;[m
[32m+[m			[32mdomi.Particles[i]->vb(2)		= -vcompress;[m
[32m+[m[32m//			domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);[m
[32m+[m		[32m}[m
[32m+[m		[32mif (domi.Particles[i]->ID == 2)[m
[32m+[m		[32m{[m
[32m+[m			[32mdomi.Particles[i]->a		= Vec3_t(0.0,0.0,0.0);[m
[32m+[m			[32mdomi.Particles[i]->v		= Vec3_t(0.0,0.0,0.0);[m
[32m+[m			[32mdomi.Particles[i]->vb		= Vec3_t(0.0,0.0,0.0);[m
[32m+[m			[32m//domi.Particles[i]->VXSPH	= Vec3_t(0.0,0.0,0.0);[m
[32m+[m		[32m}[m
[32m+[m	[32m}[m
[32m+[m[32m}[m
[32m+[m
[32m+[m
[32m+[m[32musing std::cout;[m
[32m+[m[32musing std::endl;[m
[32m+[m
[32m+[m[32mint main(int argc, char **argv) try[m
[32m+[m[32m{[m
[32m+[m[32m  SPH::Domain	dom;[m
[32m+[m
[32m+[m[32m  dom.Dimension	= 3;[m
[32m+[m[32m  dom.Nproc	= 4;[m
[32m+[m[32m  dom.Kernel_Set(Qubic_Spline);[m
[32m+[m[32m  //dom.Kernel_Set(Hyperbolic_Spline);[m
[32m+[m[32m  dom.Scheme	= 1;	//Mod Verlet[m
[32m+[m[32m  //dom.XSPH	= 0.1; //Very important[m
[32m+[m
[32m+[m[32m  double dx,h,rho,K,G,Cs,Fy;[m
[32m+[m[32m  double R,L,n;[m
[32m+[m
[32m+[m[32m  R	= 0.15;[m
[32m+[m[32m  L	= 0.56;[m
[32m+[m[32m  n	= 30.0;		//in length, radius is same distance[m
[32m+[m
[32m+[m[32m  rho	= 2700.0;[m
[32m+[m[32m  K	= 6.7549e10;[m
[32m+[m[32m  G	= 2.5902e10;[m
[32m+[m[32m  Fy	= 300.e6;[m
[32m+[m[32m  //dx	= L / (n-1);[m
[32m+[m[32m  //dx = L/(n-1);[m
[32m+[m[32m  dx = 0.015;[m
[32m+[m[32m  h	= dx*1.2; //Very important[m
[32m+[m[32m    Cs	= sqrt(K/rho);[m
[32m+[m
[32m+[m[32m    double timestep;[m
[32m+[m[32m    timestep = (0.2*h/(Cs));[m
[32m+[m
[32m+[m[32m  //timestep = 2.5e-6;[m
[32m+[m
[32m+[m[32m  cout<<"t  = "<<timestep<<endl;[m
[32m+[m[32m  cout<<"Cs = "<<Cs<<endl;[m
[32m+[m[32m  cout<<"K  = "<<K<<endl;[m
[32m+[m[32m  cout<<"