PersianSPH - Open Programming Library for Mechanical Systems

Copyright (C) 2013 Maziar Gholami Korzani and Sergio Galindo Torres

PersianSPH is a programming library for the implementation of simulation 
tools in mechanics. Its source code is mainly written in C++ with easy
to use templates for further custumization.

Potential applications include, for instance, solid and soil mechanics 
and fluid dynamics using Smoothed Particle Hydrodynamics method

For more information, please check http://korzani.wix.com/persiansph

We're going to need 4 libraries:
BLITZ++   --- VERSION 0.9!!! NOT Current 	https://github.com/luchete80/blitz-0.9-cmake.git
HDF5		
GSL											https://github.com/ampl/gsl/tags
LAPACK

VISUAL STUDIO CONFIG

CompactNSearch and CuNSearch must be separated (either they are in conflict)

1)  For blitz se debe copiar el archivo que hay en Blitz-VS.NET.zip
    al binario path%binario\blitz\ms
2) 
3) GSL: Obtained from here: 
4) The files for visual studio are located in :

Thermal - Struct Solver Steps

StartAcceleration(); ----->THIS RESETS STRAIN RATE!!!!
PrimaryComputeAcceleration();
LastComputeAcceleration(); //--->>>CALCULATES Strain Rate BUT FROM ZERO!
														//AND CALCULATES ACCELERATION FROM SIGMA (PREV. IT)

CalcConvHeat();
CalcTempInc();
CalcThermalExpStrainRate(); //Add Thermal expansion Strain Rate Term
//...... BUT SIGMA IS ALREADY CALCULATED!
Move(deltat); /// Uses ACCELERATION TO MOVE PARTICLE (calc u && v)
//CALCULATE SIGMA FOM STRAIN RATE