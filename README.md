# WeldForm

## What is it? 
WeldForm SPH is a CPU based Smooth Particle Hydrodynamics solver applied to high deformation model metal forming and processes.
GPU Version is under development here 

Has been adapted to work on both Linux Windows system (both on MinGW and MSVC compilers).

![alt text](https://github.com/luchete80/WeldForm/blob/master/compression.PNG)

Is hevaily based on: 

1) PersianSPH - Maziar Gholami Korzani and Sergio Galindo Torres (https://github.com/mghkorzani/persiansph)
2) Kirk Fraser Thesis on SPH model of FSW - https://constellation.uqac.ca/4246/1/Fraser_uqac_0862D_10345.pdf

## Features
Has been exclusively adapted to solid mechaincs, and it includes:

- Mechanic Solver
- Thermal Solver
- Coupled ThermoMechanical Solver (in progress)
- Contact formulation (in progress)
- Adaptive search only in case of plastic strain threshold (in progress)


We're going to need 4 libraries:
1 - BLITZ++   --- VERSION 0.9!!! NOT Current 	https://github.com/luchete80/blitz-0.9-cmake.git
2 - HDF5		
3 - GSL											https://github.com/ampl/gsl/tags
4 - LAPACK

GCC Config (Working on Linux and MinGW)
1 - Download precompiled libraries here: 
2 - 


VISUAL STUDIO CONFIG

CompactNSearch and CuNSearch must be separated (either they are in conflict)

1)  For blitz se debe copiar el archivo que hay en Blitz-VS.NET.zip
    al binario path%binario\blitz\ms
2) 
3) GSL: Obtained from here: 
4) The files for visual studio are located in :

