# WeldForm

WeldForm SPH is a CPU based Smooth Particle Hydrodynamics solver applied to high deformation model metal forming and processes.
Similar Projects
[GPU Version](https://github.com/luchete80/WeldFormGPU)
[WeldFormFEM](https://github.com/luchete80/WeldFormFEM) Explicit FEM all-in-one CPU/GPU Solver with Adaptive Mesh Refinement / ALE 

Tutorials:
[My channel](https://www.youtube.com/@opensourcemechanics)

Windows Binaries: https://sourceforge.net/projects/weldform/files/


Has been adapted to work on both Linux and Windows system (both on MinGW and MSVC compilers).

![alt text](https://github.com/luchete80/WeldForm/blob/master/compression_axil.png)

![alt text](https://github.com/luchete80/WeldForm/blob/master/compression.PNG)

![alt text](https://github.com/luchete80/WeldForm/blob/master/met_form.png)



Is hevaily based on: 

1) [PersianSPH](https://github.com/mghkorzani/persiansph) - Maziar Gholami Korzani and Sergio Galindo Torres
2) Kirk Fraser [Thesis](https://constellation.uqac.ca/4246/1/Fraser_uqac_0862D_10345.pdf) and [works](https://pdfs.semanticscholar.org/b09e/8c8023d56b130cc6fa5314cb66bce364df8e.pdf) on SPH model of FSW

## Features
Has been exclusively adapted to solid mechaincs, and it includes:

- Mechanic Solver
- Thermal Solver
- Coupled ThermoMechanical Solver (in progress)
- Contact formulation (in progress)
- Adaptive search only in case of plastic strain threshold (in progress)
- Parallel Reduction Sums by standard thread locking or with Nishimura table. 

## Building Instructions

We're going to need 4 libraries:
1 - BLITZ++   --- VERSION 0.9!!! (NOT Current/Last version) 	https://github.com/luchete80/blitz-0.9-cmake.git
2 - HDF5		
3 - GSL https://github.com/ampl/gsl/tags
4 - LAPACK (Not necesary by now)

For each library compile and install in PKG dir, with -DCMAKE_INSTALL_PREFIX=$PKG/libxxx

GCC Config (Working on Linux and MinGW):
- 1 - Download precompiled libraries [here](https://drive.google.com/drive/folders/16FoY47D_TQOd_0Cb_1ltLr--S6OVOaSi?usp=sharing): 
- 2 - Write or use an existing example
- 3 - Set environmental vars:
   >>> 
   set SPH=D:/Luciano/Numerico/WeldForm
   
   set PKG=D:/Luciano/Numerico/Libs
   >>>

- 4 - Create an empty dir and execute cmake (WeldFormDir)


VISUAL STUDIO CONFIG

CompactNSearch and CuNSearch must be separated (either they are in conflict)

1)  For blitz se debe copiar el archivo que hay en Blitz-VS.NET.zip
    al binario path%binario\blitz\ms
2) 
3) GSL: Obtained from here: 
4) The files for visual studio are located in :


## References
 * [PersianSPH] (https://github.com/mghkorzani/persiansph)
 * Kirk Fraser, ROBUST AND EFFICIENT MESHFREE SOLID THERMO-MECHANICS SIMULATION OF FRICTION STIR WELDING
 * Daisuke Nishiura , Hide Sakaguchi. Parallel-vector algorithms for particle simulations on shared-memory multiprocessors. Journal of Computational Physics 230 (2011) 1923–1938
