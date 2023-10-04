# Changelog

## [0.4.3] - TBR

### Input File 
 - Added Angular Velocity reading
 
## [0.4.2.1] - 20230822 - Hot Fixes

### Solver 
 - Deprecated KickDrift (was erroneously the default coupled thermal-mechanical)
   
### Input File 
 - Now "Thermal" is default Leapfrog thermal mechanical
 
## [0.4.2] - 20230818 - Bug Fixed and Improvements

### Solver 
 - CORRECTED Randles and Libersky Central Differences (Leapfrog form) solver. Now 
   accepts as Fraser Solver (fixed ts) a CFL up until 1 (default 0.6) in some cases.
   This is a MAJOR ISSUE since allows for same CFL as Fraser, but works also
   with variable step size (Fraser works only with constant time steps).
   Changed default to this solver.
   
### Input File 
 - Fixed Amplitude reading (it was an error).
 - Added amplitude function to standard velocity BC (no contact)
 
  ### Deprecated
 - KickDrift Randles & Libersky Solver (which was in Verlet Form), was a low CFL 
 
## [0.4.1] - 20230725

### Solver 
  - Now reading constant thermal properties and initial conditions. First thermal-mechanical example ready.
  - Fixed Friction Stir Welding example (greater particle density results in not so negative temp at the begining)
  - Added input reading details (to know in which vars exists an input error)
  - Randles and Libersky variable time step integrator working
### Input File 
  - Nb frequency update option on input 
  - Changed Sjohnson Cook constant order
  - Added Cylindrical grid option in input 
  ### Deprecated
 - Original Simultaneous Solver (All vars calculated at same time)

## [0.4.0] - 2023-02-03
### Core 
 - Structural , thermal and coupled solver.
 - Fraser, Leapfrog, Verlet and Differential KickDrift time integration working.
 - Bilinear, Hollomon and Johnson Cook materials.
 - Frictional Wang contact.
 - Working accumulated frictional contact Work, external forces Work and plastic dissipated Work.  
 
### Input File 
 - Wang Frictional contact algorithm working, giving same results that FEM solvers.
 - First imput file reader: working only for structural problems.
 - Input Boundary conditions of velocity and zones. 
 


