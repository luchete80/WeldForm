# WeldForm

WeldForm SPH is a CPU based Smooth Particle Hydrodynamics solver applied to high deformation model metal forming and processes.

Windows Binaries: https://sourceforge.net/projects/weldform/files/
GPU Version is under development [here](https://github.com/luchete80/WeldFormGPU)

Has been adapted to work on both Linux and Windows system (both on MinGW and MSVC compilers).

![alt text](https://github.com/luchete80/WeldForm/blob/master/compression.PNG)

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

## Building Instructions

We're going to need 4 libraries:
1 - BLITZ++   --- VERSION 0.9!!! (NOT Current/Last version) 	https://github.com/luchete80/blitz-0.9-cmake.git
2 - HDF5		
3 - GSL											https://github.com/ampl/gsl/tags
4 - LAPACK

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

